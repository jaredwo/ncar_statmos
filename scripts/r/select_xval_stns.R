# Script to build a longform dataframe of precipitation observations and covariates
# for input to the neural network interpolation method

library(rgdal)
library(reshape2)
# Set timezone for environment to UTC
Sys.setenv(TZ="UTC")

# Constant for path to local ncar statmos code location
PATH_CODE <- '/storage/home/jwo118/workspace/ncar_statmos'
# Constant for path to local netcdf data file
FPATH_DATA <- '/storage/home/jwo118/scratch/ncar_statmos/data/prcp_ghcnd_midatlantic_19480101_20151231.nc'
START_TIME = '2011-01-01'
END_TIME <- '2015-12-31'
FPATH_DF <- '/storage/home/jwo118/scratch/ncar_statmos/data/prcp_stns_cheasapeake_20110101_20151231.rdata'
load(FPATH_DF)
stnids <- sort(unique(as.character(df_prcp$station_id)))

# Source in utility functions and packages for loading station data
source(file.path(PATH_CODE, 'interp_r', 'util_data.R'))

# Open netcdf observation file
ds <- h5file(FPATH_DATA, mode='r')

# Load station metadata as a SpatialPointsDataFrame
stns <- obsnc_stns_spdf(ds)
stns <- stns[stnids,]

# Load time variable data and create vector of POSIXct times
# This is needed when loading actual observations
times <- obsnc_time(ds)

# Load precipitation time series for stations
xts_prcp <- obsnc_obs(ds, 'prcp', stns, times, start_end=c(START_TIME, END_TIME), stnids=stns$station_id)

n_xval <- round(length(stnids)*.1)
mask_por <- (colSums(is.finite(xts_prcp)) / nrow(xts_prcp)) >= .9
stnids_xval <- stns$station_id[mask_por]
stnids_xval <- sample(stnids_xval,n_xval)

save(stnids_xval,file='/storage/home/jwo118/scratch/ncar_statmos/data/stnids_xval.rdata')
