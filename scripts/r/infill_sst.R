# Multiprocessing script for infilling missing data in the CoRTAD data product
# using universal kriging.

# Load required libraries
library(raster)
library(ncdf4)
library(gstat)
library(sp)
# Libaries for multiprocessing
library(doParallel)
library(foreach)

# Constants
CRS_WGS84 <- CRS('+proj=longlat +datum=WGS84')
FPATH_SST_NC <- '[Path to the input CoRTAD netcdf file]'
PATH_OUT <- '[Path to ouput directory to place infilled GeoTIFF raster files]'
FPATH_LOGFILE <- '[Path to log file where progress will be written]'
NPROCS <- 10 #number of processors

# Get land/sea mask
dslndm <- raster(FPATH_SST_NC,varname='land')
proj4string(ds) <- CRS_WGS84
proj4string(dslndm) <- CRS_WGS84
date_names <- names(ds)
mask_sea <- dslndm[] == 0
rm(dslndm)
gc()

# Cluster setups

fpath_log <- FPATH_LOGFILE
writeLines(c(""),fpath_log)
cl <- makeCluster(NPROCS)
registerDoParallel(cl)

# Multiprocessing loop
results <- foreach(a_date=date_names) %dopar% { 
	
	# Need to reload libraries in cluster child process
	library(raster)
	library(ncdf4)
	library(gstat)
	library(sp)
	
	log_file <-file(fpath_log,open='a')
	writeLines(sprintf("Started processing %s", as.character(a_date)),log_file)
	flush(log_file)
	
	ds = brick(FPATH_SST_NC)
	proj4string(ds) <- CRS_WGS84
	
	sgdf <- as(ds[[a_date]],'SpatialGridDataFrame')
	spdf <- as(sgdf,'SpatialPointsDataFrame')
	spdf$x <- spdf@coords[,1]
	spdf$y <- spdf@coords[,2]
	
	interp_grid <- ds[[a_date]]
	lon_grid <- init(interp_grid,'x')
	lat_grid <- init(interp_grid,'y')
	
	mask_interp <- is.na(interp_grid[]) & mask_sea
	lon_grid[!mask_interp] <- NA
	lat_grid[!mask_interp] <- NA
	interp_stack <- stack(list(x=lon_grid,y=lat_grid))
	spgdf_interp <- as(interp_stack,'SpatialGridDataFrame')
	
	spdf_sample <- spdf[sampleInt(nrow(spdf),25000),]
	evgm <- variogram(as.formula(sprintf("%s~x+y",a_date)),spdf_sample,cutoff=6000)
	a_vgm <- fit.variogram(evgm,vgm(c('Gau','Sph','Exp')), fit.kappa=TRUE)
	
	sst_interp <- krige(as.formula(sprintf("%s~x+y",a_date)),spdf_sample,spgdf_interp,model=a_vgm,nmax=30)
	sst_interp <- raster(sst_interp['var1.pred'])
	interp_grid[mask_interp] <- sst_interp[mask_interp]
	
	fname_out <- sprintf('sst_%s.tif',a_date)
	writeRaster(interp_grid,file.path(PATH_OUT,fname_out))
	
	writeLines(sprintf("Finished processing %s", as.character(a_date)),log_file)
	flush(log_file)
}
stopCluster(cl)
