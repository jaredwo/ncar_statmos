# Example of using Random Classification Forests and Quantile Random Regression
# Forests to interpolate precipitation occurence and amount on a single day.
###############################################################################

Sys.setenv(TZ="UTC")

# Constants for script
PATH_CODE <- '[PLACE PATH TO CODE DIRECTORY HERE]'
FPATH_STNDATA <- '[PLACE FILE PATH TO NETCDF FILE HERE]'
DATE_INTERP <- as.Date('1980-05-14')
NSTNS_XVAL <- 150
NTREE_RF <- 500

# Source in utility functions and packages for loading station data
source(file.path(PATH_CODE, 'interp_r', 'util_data.R'))
# Source in utility functions and packages for interpolation
source(file.path(PATH_CODE, 'interp_r', 'util_interp.R'))

################################################################################
# Load station data
################################################################################

# Open netcdf observation file
ds <- h5file(FPATH_STNDATA, mode='r')

# Load station metadata as SpatialPointsDataFrame
stns <- obsnc_stns_spdf(ds)
# Add coordinates as actual variables on the DataFrame
stns$latitude1 <- stns@coords[,'latitude']
stns$longitude1 <- stns@coords[,'longitude']

# Load dates for the time series
times <- obsnc_time(ds)
# Load observations for all stations (stnids=NULL) for the day to interpolate.
# Observations are returned as an xts object
obs <- obsnc_obs(ds, 'prcp', stns, times, start_end=c(as.character(DATE_INTERP),
				 as.character(DATE_INTERP)), stnids=NULL)
 
# Add observations for the day as a prcp variable to stns
stns$prcp <- as.numeric(obs)
 
# Remove duplicate stations
not_dup_mask <- !as.logical(ds['dup'][])

stns <- stns[not_dup_mask,]

# Close netcdf observation file
h5close(ds)

# Drop stations that do not have an observation for the day
stns <- stns[is.finite(stns$prcp),]

# Add wet variable
stns$wet <- stns$prcp > 0

# Select and drop n random cross validation stations to be used for validation
# To make sure random validation stations are spread across the domain, cluster
# by lon/lat location and then select a station from each cluster
clusters <- kmeans(as.data.frame(stns)[,c('longitude1','latitude1')],
				   centers=NSTNS_XVAL)$cluster
stnids_xval <- as.character(aggregate(stns$station_id, list(cluster=clusters),
				function(x) x[sample(seq_along(x),1)])$x)
stnids_xval <- sort(stnids_xval)
stns_xval <-  stns[stnids_xval,]
stns <- stns[!(stns$station_id %in% stnids_xval),]

################################################################################
# Interpolate precipitation occurrence
################################################################################
xval_rf <- kfold_xval(stns,'wet',stns,fold=5,pred_func=func_rf_wet,a_seed=76, ntree=NTREE_RF)
xval_gam <- kfold_xval(stns,'wet',stns,fold=5,pred_func=func_gam_wet,a_seed=76)

aucs <- tryCatch(c(auc(xval_rf$wet,xval_rf$wet_p),auc(xval_gam$wet,xval_gam$wet_p)), error=function(e) c(.5,.5))
wetts <- c(optimal_threshold(xval_rf$wet,xval_rf$wet_p),optimal_threshold(xval_gam$wet,xval_gam$wet_p))
wgts <- aucs/sum(aucs)

pred_rf <- func_rf_wet(stns, newdata=stns_xval, ntree=NTREE_RF)
pred_gam <- func_gam_wet(stns,newdata=stns_xval)
pred_combined <- (pred_rf*wgts[1]) + (pred_gam*wgts[2])
pred_combined <- pred_combined > weighted.mean(wetts,wgts)
pred_wet <- pred_combined

################################################################################
# Precipitation Amount
################################################################################

stns_wet <- stns[stns$wet,]
stns_xval_wet <- stns_xval[pred_wet,]

pred_amt <- func_rfq_amt(stns_wet,stns_xval_wet,NTREE_RF)

#Add final precipitation predictions as variable to stns_xval
pred_fnl <- rep(0,nrow(stns_xval))
pred_fnl[pred_wet] <- pred_amt
stns_xval$prcp_p <- pred_fnl