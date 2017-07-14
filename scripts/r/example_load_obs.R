# Example script for loading station precipitation observations from a netcdf file.
# Requires the following packages:
#
# h5: provides utilities for working with netcdf4/hdf5 files
# RNetCDF: provides utilities for working with netcdf files
# ncdf4: provides utilities for working with netcdf files
# xts: provides temporal data structures
# sp: provides spatial data structures
# spacetime: provides spatiotemporal data structures
# reshape2: provides utilities for wrangling data
# stringr: provides utilities for string manipulation

# Constant for path to local ncar statmos code location
PATH_CODE <- '[PLACE PATH TO CODE DIRECTORY HERE]'
# Constant for path to local netcdf data file
FPATH_DATA <- '[PLACE FILE PATH TO NETCDF FILE HERE]'

# Set timezone for environment to UTC
Sys.setenv(TZ="UTC")

# Source in utility functions and packages for loading station data
source(file.path(PATH_CODE, 'interp_r', 'util_data.R'))

# Open netcdf observation file
ds <- h5file(FPATH_DATA, mode='r')

# Load station metadata as a SpatialPointsDataFrame
stns <- obsnc_stns_spdf(ds)

# Plot the station points
plot(stns)

# Load time variable data and create vector of POSIXct times
# This is needed when loading actual observations
times <- obsnc_time(ds)

# Load precipitation time series for Richmond airport (station id: GHCND_USW00013740)
# Observations are returned as an xts object
ts_prcp <- obsnc_obs(ds, 'prcp', stns, times, stnids=c('GHCND_USW00013740'))

# Plot Richmond time series
plot(ts_prcp)

# Load observations from all stations for May 1980
# Observations are returned as an xts object with a column for each station
# This is considered a "spacewide" data structure
prcp <- obsnc_obs(ds, 'prcp', stns, times, start_end=c('1980-05-01','1980-05-31'))

# Convert spacewide xts object to a spatio-temporal data structure from the
# spacetime package
prcp_stfdf <- coerce_spacewide_xts_to_STFDF(prcp, stns, 'prcp')

# Plot map of observations for May 14, 1980
spplot(prcp_stfdf[,'1980-05-14','prcp'])

# Plot multiple maps of observations from May 14-21
stplot(prcp_stfdf[,'1980-05-14/1980-05-21','prcp'])
