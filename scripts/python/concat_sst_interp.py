from interp.raster import open_as_da #
import glob
import os
import pandas as pd
import xarray as xr

PATH_INFILLED_RATER = '[Path to infilled raster GeoTIFFs]'
FPATH_OUT = '[File path for final output netcdf file]'

if __name__ == '__main__':
    
    fpaths = sorted(glob.glob(os.path.join(PATH_INFILLED_RATER, 'sst_*.tif' )))
    da_all = []
    dates = []
    
    for fpath in fpaths:
        
        da = open_as_da(fpath)
        da_all.append(da)
        dates.append(pd.Timestamp(os.path.basename(fpath)[5:15]))
    
    itime = pd.DatetimeIndex(dates,name='time')
    ds = xr.concat(da_all,dim=itime)
    ds.to_dataset(name='sst').to_netcdf(FPATH_OUT)