'''
# Example script for loading station precipitation observations from a netcdf file.
# Requires the follow packages:
# xarray (http://xarray.pydata.org)
# numpy (http://www.numpy.org)
# pandas (http://pandas.pydata.org)
# netCDF4 (http://www.numpy.org)
# matplotlib (https://matplotlib.org)
# cartopy (http://scitools.org.uk/cartopy)
'''

# Import required libraries
import sys
import xarray as xr

# Constant for path to local ncar statmos code location
PATH_CODE = '[PLACE PATH TO CODE DIRECTORY HERE]'

# Constant for path to local netcdf data file
FPATH_DATA = '[PLACE FILE PATH TO NETCDF FILE HERE]'

# Append code path to python path and import map plotting utility function
sys.path.append(PATH_CODE)
from interp.util_data import map_pts

if __name__ == '__main__':
    
    # Open netcdf file of station precipitation observations
    ds = xr.open_dataset(FPATH_DATA)
    
    # Load station metadata as pandas DataFrame
    stns = ds[['station_name','elevation','longitude','latitude']].to_dataframe()
    
    # Plot the station points
    map_pts(stns.longitude.values, stns.latitude.values)
    
    # Plot precipitation time series for a single station
    # Richmond airport (station id: GHCND_USW00013740)
    ds.prcp.loc[:,'GHCND_USW00013740'].plot()
    
    # Load observation totals for May 21, 1980 as pandas DataFrame
    obs = ds.prcp.loc['1980-05-21'].to_dataframe().loc[:,'prcp']
    
    # Join with stns DataFrame
    obs = stns.join(obs)
    
    # Drop stns with missing observations
    obs = obs.dropna()
    
    # Plot observation totals
    map_pts(obs.longitude.values, obs.latitude.values, colors=obs.prcp.values)
    