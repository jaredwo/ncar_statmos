# Functions for loading/writing station data from/to netCDF observation files
# along with general data utility functions.
###############################################################################

library(h5)
library(RNetCDF)
library(ncdf4)
library(xts)
library(sp)
library(reshape2)
library(stringr)
library(spacetime)

#' Load station metadata and create station SpatialPointsDataFrame
#' 
#' @param ds H5File object pointing to netCDF observation file
#' @param vnames variables names of station metadata to load
#' @param vname_stnid variable name from vnames representing unique station identifier
#' @param vnames_xy variable names from vnames representing x,y position
#' @return SpatialPointsDataFrame
obsnc_stns_spdf <- function(ds, vnames=c('station_id','station_name',
                                         'longitude','latitude', 'elevation'),
                            vname_stnid = 'station_id',
                            vnames_xy = c('longitude','latitude')) {
  
  stns <- list()
  for (vname in vnames) {
    
    vals <- ds[vname][]
    
    # Create vector of strings if 2D character
    if (is.character(vals) & length(dim(vals)) == 2) {
      vals <- apply(vals,1,paste0,collapse="")
    }
    
    stns[[vname]] <- vals
  }
  stns[['obs_index']] <- 1:length(stns[[1]])
  
  # Turn station dataframe into spatial dataframe
  stns <- as.data.frame(stns, stringsAsFactors = FALSE)
  row.names(stns) <- stns[,vname_stnid]
  coordinates(stns) <- as.formula(paste0("~",vnames_xy[1],"+",vnames_xy[2]))
  proj4string(stns)=CRS("+proj=longlat +datum=WGS84")
  
  return(stns)
  
}

#' Load time variable data and create vector of POSIXct times
#' 
#' @param ds H5File object pointing to netCDF observation file
#' @param vname_time variable name representing time
#' @return POSIXct vector of times
obsnc_time <- function(ds, vname_time='time') {
  
  times <- ds[vname_time][]
  times <- utcal.nc(h5attr(ds[vname_time],'units'), times, type="c")
  return(times)
 
}

#' Load station observations
#' 
#' @param ds H5File object pointing to netCDF observation file
#' @param vname variable name for which to load observations
#' @param stns SpatialPointsDataFrame of station metadata
#' @param times POSIXct vector of times for the netCDF observation file
#' @param start_end character vector of size 2 with start and end dates
#' @param stnids station ids for which to load observations
#' @return xts object of observations
obsnc_obs <- function(ds, vname, stns, times, start_end=NULL, stnids=NULL) {
  
  if (is.null(start_end) & is.null(stnids))
  {
    
    stnids <- stns$station_id  
    obs <- ds[vname][]
    
  } else if (is.null(start_end) & (!is.null(stnids))) {
    
    i <- stns@data[stnids, 'obs_index']
    
    if (is.unsorted( i, strictly = TRUE))
      stop('Station ids must be in obs_index order')
    
    obs <- ds[vname][,i]
    
  } else if ((!is.null(start_end)) & is.null(stnids)) {
    
    stnids <- stns$station_id  
    i <- which(times >= start_end[1] & times <= start_end[2])
    times <- times[i]
    obs <- ds[vname][i,]
    
  } else {
    i_stns <- stns@data[stnids, 'obs_index']
    
    if (is.unsorted( i_stns, strictly = TRUE))
      stop('Station ids must be in obs_index order')
    
    i_time <- which(times >= start_end[1] & times <= start_end[2])
    times <- times[i_time]
    obs <- ds[vname][i_time,i_stns]
  }
  
  colnames(obs) <- stnids
  
  if ("missing_value" %in% list.attributes(ds[vname])) {
    obs[obs==h5attr(ds[vname], 'missing_value')] = NA
  }
  
  obs <- xts(obs, times)
  
  return(obs)
  
}

#' Build a boolean mask for stations that have long enough period-of-record
#' 
#' @param ds H5File object pointing to netCDF observation file
#' @param vname variable name for which to build the mask
#' @param start_date The start date for period-of-record time period 
#' @param end_date The end date for period-of-record time period
#' @param nyrs The minimum period of record in years. The function tests
#' whether a station has at least nyrs years of data in each month. 
#' @return boolean mask for stations 
obsnc_pormask <- function(ds, vname, start_date, end_date, nyrs) {
  
  vname <- paste0('obs_cnt_prcp_', format.Date(start_date,'%Y%m%d'),'_',
                  format.Date(end_date,'%Y%m%d'))
  
  obs_cnts <- ds[vname][]
  
  mths <- as.POSIXlt(seq(as.Date('2015-01-01'),
                     as.Date('2015-12-31'), by='day'))$mon + 1
  
  days_in_month <- sapply(1:12,FUN=function(x){sum(mths==x)})
  
  nmin <- days_in_month*nyrs
  
  mask_por <- t(sapply(1:12,function(x){obs_cnts[x,] > nmin[x]}))
  
  mask_por <- colSums(mask_por) == 12
  
  return(mask_por)
  
}

#' Initialize a new netCDF observation file
#' @param fpath file path for the new netCDF observation file 
#' @param stns SpatialPointsDataFrame with station metadata
#' @param times POSIXct vector of observations times
#' @param vnames vector observation variable names (e.g.--tmin, tmax, prcp)
#' @return H5File object pointing to new netCDF observation file
obsnc_init <- function(fpath, stns, times, vnames) {
 
  # Create station_id dimension
  dim_stnid  <- ncdim_def('station_id', '', seq(length(stns$station_id)),
                            create_dimvar=FALSE)
  
  # Create time dimension
  times <- as.POSIXct(times)
  tunits <- sprintf("days since %s", strftime(min(times),"%Y-%m-%d"))
  dim_time <- ncdim_def('time', tunits, utinvcal.nc(tunits, times))
  
  # Create station metadata variables
  vars_meta <- list()
  vals_meta <- list()
  i <- 1
  
  # First loop through data variables
  for (vname in colnames(stns@data)) {
    
    a_vals_meta <- stns@data[,vname]
    
    if (class(a_vals_meta) == 'character') {
      
      # If character datatype, values must be stored in 2d variable
      char_len <- max(nchar(a_vals_meta))
                  
      # Create string dimension with size char_len
      dim_string  <- ncdim_def(paste0("string_",vname), '', seq(char_len),
                               create_dimvar=FALSE)
      # Create 2-d variable
      a_var <- ncvar_def(vname, '', list(dim_string, dim_stnid), prec='char')
      
    } else {
      
      a_var <- ncvar_def(vname, '', list(dim_stnid), prec='double')
      
    }
    
    vars_meta[[i]] <- a_var
    vals_meta[[i]] <- a_vals_meta
    i = i + 1
    
  }
  
  # Add coordinate variables to the metadata
  for (vname in colnames(coordinates(stns))) {
  
    a_vals <- as.numeric(coordinates(stns)[, vname])
    a_var <- ncvar_def(vname, '', list(dim_stnid), prec='double')
    
    vars_meta[[i]] <- a_var
    vals_meta[[i]] <- a_vals
    i = i + 1
  
  }
  
  # Create main 2D observation variables
  vars_main = list()
  i = 1
  
  for (vname in vnames) {
    a_var <- ncvar_def(vname, '', list(dim_stnid, dim_time), prec='double',
                       compression=4, chunksizes=c(1,length(times)), missval=NA)
    vars_main[[i]] <- a_var
    i = i + 1
  }
  
  # Create netcdf file and add station metadata
  ds <- nc_create(fpath, append(vars_meta,vars_main), force_v4=T)
  
  for (i in seq(length(vals_meta))) {
    
    a_var <- vars_meta[[i]]
    a_vals_meta <- vals_meta[[i]]
    
    ncvar_put(ds, a_var, a_vals_meta)    
    
  }
  
  nc_sync(ds)
  nc_close(ds)
  
  # Return h5file object pointing to the new netcdf file
  return(h5file(fpath))
  
}

#' Convert a tidy/longform dataframe into a spacewide xts object
#' @param df_tidy 
#' @param value.var 
#' @param time_col 
#' @param id_col 
#' @return 
coerce_tidy_df_to_spacewide_xts <- function(df_tidy, value.var, time_col='time',
                                            id_col='station_id') {
  
  xts_sw <- dcast(df_tidy, as.formula(paste0(time_col,"~",id_col)), value.var=value.var)
  row.names(xts_sw) <- xts_sw[[time_col]]
  xts_sw <- as.xts(xts_sw[,colnames(xts_sw) != time_col])   
  return(xts_sw)
   
}

#' Convert a spacewide xts object to a spacetime STFDF object
#' @param xts_sw 
#' @param spdf_locs 
#' @param varname 
#' @return 
coerce_spacewide_xts_to_STFDF <- function(xts_sw, spdf_locs, varname) {
  
  # Create spacetime STFDF object.
  # Data needs to be a single column data.frame with station id varying the fastest
  # By default time varies the fastest with as.vector call, so transpose first
  
  adf <- data.frame(as.vector(t(xts_sw)))
  colnames(adf) <- varname  
  a_stfdf <- STFDF(spdf_locs, index(xts_sw), data=adf)
  
  return(a_stfdf)
}