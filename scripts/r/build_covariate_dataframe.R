# Script to build a longform dataframe of precipitation observations and covariates
# for input to the neural network interpolation method

library(rgdal)
library(reshape2)

# Constant for path to local ncar statmos code location
PATH_CODE <- '/storage/home/jwo118/workspace/ncar_statmos'
# Constant for path to local netcdf data file
FPATH_DATA <- '/storage/home/jwo118/scratch/ncar_statmos/data/prcp_ghcnd_midatlantic_19480101_20151231.nc'
# Constant for path to Chesapeake Bay watershed shapefile
PATH_CBAY <- '/storage/home/jwo118/scratch/topowx-prcp/vector/'
# Start/end date constants
START_TIME = '2006-01-01'
END_TIME <- '2015-12-31'
# Constant for output rdata file
FPATH_OUT <- '/storage/home/jwo118/scratch/ncar_statmos/data/prcp_stns_cheasapeake_20060101_20151231.rdata'

# Set timezone for environment to UTC
Sys.setenv(TZ="UTC")

# Source in utility functions and packages for loading station data
source(file.path(PATH_CODE, 'interp_r', 'util_data.R'))

# Open netcdf observation file
ds <- h5file(FPATH_DATA, mode='r')

# Load station metadata as a SpatialPointsDataFrame
stns <- obsnc_stns_spdf(ds)

# Subset stations to just the Chesapeake Bay watershed
# Read Chesapeake Bay watershed polygon
poly_cbay <- as(readOGR(PATH_CBAY, 'chesapeake_basin_nad83'), 'SpatialPolygons')
# Make coordinate system same as stations (NAD83 vs. WGS84 differences shouldn't matter)
proj4string(poly_cbay) <- proj4string(stns)
mask_cbay <- is.finite(over(stns,poly_cbay))
stns <- stns[mask_cbay,]

# Load time variable data and create vector of POSIXct times
# This is needed when loading actual observations
times <- obsnc_time(ds)

# Load precipitation time series for stations
xts_prcp <- obsnc_obs(ds, 'prcp', stns, times, start_end=c(START_TIME, END_TIME), stnids=stns$station_id)
yday <- .indexyday(ts_prcp) + 1
df_prcp <- as.data.frame(xts_prcp)
df_prcp$time <- index(xts_prcp)
df_prcp$yday <- .indexyday(xts_prcp) + 1

# Melt into long form dataframe and drop NAs
df_prcp <- melt(df_prcp,id.vars=c('time','yday'), na.rm=TRUE, value.name='prcp', stringsAsFactors=FALSE)
names(df_prcp)[names(df_prcp) == 'variable'] <- 'station_id'

# Join with stns dataframe to add lon, lat, and elev to each observation
df_stns <- as.data.frame(stns,stringsAsFactors=FALSE)[,c('station_id','longitude','latitude','elevation')]
df_prcp <- merge(df_prcp, df_stns, by='station_id', all.x=TRUE)

# Add wet variable
df_prcp$wet <- df_prcp$prcp > 0

# Add cos/sin seasonality variables
df_prcp$yday_sin <- sin((2*pi*df_prcp$yday)/365.25)
df_prcp$yday_cos <- cos((2*pi*df_prcp$yday)/365.25)

#Output dataframe as rdata file
save(df_prcp,file=FPATH_OUT)