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
START_TIME = '2011-01-01'
END_TIME <- '2015-12-31'
# Constant for output rdata file
FPATH_OUT <- '/storage/home/jwo118/scratch/ncar_statmos/data/prcp_30stns_cheasapeake_20110101_20151231.rdata'
N_NN <- 2

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

# Only keep 30 stations with longest por
#stnids_keep <- colnames(xts_prcp)[order(colSums(is.finite(xts_prcp)), decreasing=TRUE)[1:30]]
#xts_prcp <- xts_prcp[,stnids_keep]

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

# Get rid of stations that didn't have any observations for period specified
mask_has_obs <- stns$station_id %in% unique(as.character(df_prcp$station_id))
stns <- stns[mask_has_obs,]
df_stns <- df_stns[mask_has_obs,]

# Add N nearest neighbor covariates
for (i in seq(N_NN)) {
	df_prcp[sprintf('dlon_nn%.2d', i)] <- rep(NA, nrow(df_prcp))
	df_prcp[sprintf('dlat_nn%.2d', i)] <- rep(NA, nrow(df_prcp))
	df_prcp[sprintf('delev_nn%.2d', i)] <- rep(NA, nrow(df_prcp))
	df_prcp[sprintf('prcp_nn%.2d', i)] <- rep(NA, nrow(df_prcp))
}

# Distance Matrix
d <- spDists(stns,longlat=TRUE)
row.names(d) <- stns$station_id
colnames(d) <- stns$station_id

for (i in seq(nrow(stns))) {
	
	stnid <- stns$station_id[i]
	print(sprintf("%s: station %d of %d",stnid,i,nrow(stns)))
	
	d_i <- d[i,-i]
	ngh_ids <- names(d_i[order(d_i)])
	
	j_rows <- which(df_prcp$station_id==stnid)
	j_dates <- df_prcp[j_rows,'time']
	
	ngh_avail <- matrix(TRUE,nrow=length(j_rows),ncol=length(ngh_ids),dimnames=list(NULL,ngh_ids))
	
	for (n in seq(1:N_NN)) {

		for (ngh_id in ngh_ids) {
			
			prcp_ngh <- as.numeric(xts_prcp[j_dates, ngh_id])
			mask_use <- is.finite(prcp_ngh) & ngh_avail[,ngh_id] & (!is.finite(df_prcp[j_rows,sprintf('prcp_nn%.2d', n)]))
			
			if (any(mask_use)) {
				
				df_prcp[j_rows[mask_use],sprintf('dlon_nn%.2d', n)] <- df_stns[stnid,'longitude'] - df_stns[ngh_id,'longitude']
				df_prcp[j_rows[mask_use],sprintf('dlat_nn%.2d', n)] <- df_stns[stnid,'latitude'] - df_stns[ngh_id,'latitude']
				df_prcp[j_rows[mask_use],sprintf('delev_nn%.2d', n)] <- df_stns[stnid,'elevation'] - df_stns[ngh_id,'elevation']
				df_prcp[j_rows[mask_use],sprintf('prcp_nn%.2d', n)] <- prcp_ngh[mask_use]
				
				ngh_avail[mask_use,ngh_id] <- FALSE
				
			}
			
			if (all(is.finite(df_prcp[j_rows,sprintf('prcp_nn%.2d', n)]))) break
			
			
		}
		
	}

}

#Output dataframe as rdata file
save(df_prcp,file=FPATH_OUT)

# TODO: Seasonal cycle of precipitation
xts_dmean <- xts(rowMeans(xts_prcp, na.rm=TRUE),index(xts_prcp))
iyday <- .indexyday(xts_dmean) + 1

yd_mean <- rep(NA, 365)
for (i in seq(1,365)) {

	xmin <- i-45
	xmax <- i+45
	
	if (xmin < 1) {
		
		mask_xmin <- iyday >= (365 + xmin)
		mask_xmax <- iyday <= xmax
		yd_mean[i] <- mean(xts_dmean[mask_xmin | mask_xmax])
		
	} else if (xmax > 365) {
		
		mask_xmin <- iyday >= xmin
		mask_xmax <- iyday <= (xmax-365)
		yd_mean[i] <- mean(xts_dmean[mask_xmin | mask_xmax])
		
	} else {
		
		yd_mean[i] <- mean(xts_dmean[(iyday >= xmin & iyday <= xmax)])
	}
	
}
