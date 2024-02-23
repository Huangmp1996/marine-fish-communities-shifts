library(tidyverse)
# binding linear fit of maxLengthTL, Troph and tempPref
setwd('/Share/user/picapica/marine_global_change/world_ocean_climate/')
maxLengthTL_fit <- read_csv('maxLengthTL_fit.csv') %>% filter(term != '(Intercept)')
Troph_fit <- read_csv('Troph_fit.csv') %>% filter(term != '(Intercept)')
tempPref_fit <- read_csv('tempPref_fit.csv') %>% filter(term != '(Intercept)')
colnames(tempPref_fit) <- c('grid', paste0(colnames(tempPref_fit)[2:10], '.tempPref'))
sub_fit <- full_join(maxLengthTL_fit, Troph_fit, suffix = c('.maxLengthTL', '.Troph'), by = 'grid')
total_fit <- full_join(sub_fit, tempPref_fit, by = 'grid')
write.csv(total_fit, 'all_fit.csv', row.names = F)

# extract fishing effort values of grids
library(raster)
library(rasterVis)
library(exactextractr)
# csv to rasterlayer
fishing_effort <- read_csv('fishing_effort.csv') %>% 
    separate(cell_id, into = c('lat', 'lon'), sep = '_', remove = F) %>% 
    mutate(across(c(lon, lat), as.numeric))
fishing_effort_raster <- rasterFromXYZ(fishing_effort[, c('lon', 'lat', 'hours_sum')], res = c(0.1, 0.1))
rm(fishing_effort)
# saveRDS(fishing_effort_raster, 'fishing_effort_raster.rds')
fishing_effort_raster <- readRDS('fishing_effort_raster.rds')

total_fit <- read_csv('all_fit.csv')
library(dggridR)
library(sf)
dggs <- dgconstruct(res = 12, metric = FALSE)
grid <- dgcellstogrid(dggs, unique(total_fit$grid))  # shapefile grids for all fit 
fishing_effort = exact_extract(fishing_effort_raster, grid, fun = 'mean')
fishing_effort_table <- cbind(fishing_effort, grid$seqnum) %>% as.data.frame()
colnames(fishing_effort_table) <- c('fishing_effort', 'grid')

# library(ncdf4)
# ncfile <- 'merged_Hadley_NOAA//MODEL.SST.HAD187001-198110.OI198111-202206.nc'
# nc <- ncdf4::nc_open(ncfile)
# names(nc$var)
# ncatt_get(nc,'SST', "_FillValue")
# SST <- ncvar_get(nc, varid = 'SST')
# dim(SST)
# nc_close(nc)

# sst <- raster::stack(ncfile, varname = 'SST')
# empty_r <- raster(nrow = 180, ncol = 360)
# sst_stack <- stack(empty_r)
# year_cut <- seq(1, 1824, by = 12)
# for (i in year_cut){
#     year_mean <- calc(sst[[i:(i + 11)]], mean)
#     sst_stack <- addLayer(sst_stack, year_mean)
# }
# # saveRDS(sst_stack, 'sst_stack.rds')

# trend <- cellStats(sst_stack, stat='mean')
# names(trend) <- 1870:2021
# qplot(x = 1870:2021, y = trend, geom = 'line')

# extract sst values of grids between 1950-2021
sst_stack <- readRDS('sst_stack.rds')
sst_grid <- grid # sst table between 1950-2021
for (i in 81:152){
    # i = 1
    r <- rotate(sst_stack[[i]])
    sst_value <- exact_extract(r, sst_grid, fun = 'mean', progress = FALSE)
    sst_grid[ , ncol(sst_grid) + 1] <- sst_value                  # Append new column
    colnames(sst_grid)[ncol(sst_grid)] <- paste0("year_", i + 1869)  # Rename column name
}

# linear trend of sst during survey years
sst_fitting <- function(traits) {
    # traits = 'tempPref'
    sst_trend <- data.frame() # sst fitted trends
    sub_fit <- dplyr::select(total_fit, ends_with(traits), grid) %>% distinct()
    for (i in 1:nrow(sub_fit)){
        # i = 56
        year_index <- with(sub_fit, get(paste0('min_year.', traits))[i]:get(paste0('max_year.', traits))[i]) # year range of fitted grids
        sst_table <- sst_grid[which(sst_grid$seqnum == sub_fit$grid[i]), paste0('year_', year_index)] %>%
            st_set_geometry(NULL) %>% t() %>% as.data.frame()
        sst_table[ , 2] <- year_index
        colnames(sst_table) <- c('sst', 'year')
        sst_fit <- broom::tidy(lm(sst ~ year, data = sst_table))
        sst_fit$dataset_name <- with(sub_fit, rep(get(paste0('dataset_name.', traits))[i], 2))
        sst_fit$grid <- rep(sub_fit$grid[i], 2)
        sst_fit <- mutate(sst_fit, dataset_grid = paste(dataset_name, grid, sep = '+')) %>% 
            filter(term != '(Intercept)') %>% 
            dplyr::select(estimate, p.value, dataset_grid)
        colnames(sst_fit) <- paste0(colnames(sst_fit), '.', traits)
        # temp baseline
        sst_baseline <- sst_grid[which(sst_grid$seqnum == sub_fit$grid[i]), paste0('year_', c((year_index[1]-10):year_index[1]))] %>%
            st_set_geometry(NULL) %>% t() %>% as.vector() %>% mean()
        sst_fit$sst_baseline <- sst_baseline
        sst_trend <- bind_rows(sst_trend, sst_fit)
    }
    return(sst_trend)
}

sst_trend_tempPref <- sst_fitting(traits = 'tempPref')
sst_trend_maxLengthTL <- sst_fitting(traits = 'maxLengthTL')
sst_trend_Troph <- sst_fitting(traits = 'Troph')
colnames(sst_trend_tempPref) <- c("sst_estimate.tempPref", "sst_p.value.tempPref", "dataset_grid.tempPref", "sst_baseline.tempPref")
colnames(sst_trend_maxLengthTL) <- c("sst_estimate.maxLengthTL", "sst_p.value.maxLengthTL", "dataset_grid.maxLengthTL", "sst_baseline.maxLengthTL")
colnames(sst_trend_Troph) <- c("sst_estimate.Troph", "sst_p.value.Troph", "dataset_grid.Troph", "sst_baseline.Troph")

total_fit_env <- left_join(total_fit, fishing_effort_table, by = 'grid') %>% 
    left_join(sst_trend_tempPref, by = 'dataset_grid.tempPref') %>% 
    left_join(sst_trend_maxLengthTL, by = 'dataset_grid.maxLengthTL') %>% 
    left_join(sst_trend_Troph, by = 'dataset_grid.Troph')
write.csv(total_fit_env, 'all_fit_env.csv', row.names = F)
