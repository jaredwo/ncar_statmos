'''
Utilities for working with raster datasets
'''

from rasterio.warp import reproject, Resampling
import gdalconst
import netCDF4 as nc
import numpy as np
import osgeo.gdal as gdal
import osgeo.osr as osr
import pyproj
import rasterio as rio
import rasterio.mask as rio_mask
import xarray as xr

__all__ = ['RasterDataset', 'init_out_nc', 'to_inmemory_rio', 'mask_and_crop',
           'reproject_da', 'to_raster_file', 'open_as_da']

_PROJECTION_GEO_WGS84 = 4326  # EPSG Code
_PROJECTION_GEO_NAD83 = 4269  # EPSG Code
_PROJECTION_GEO_WGS72 = 4322  # EPSG Code

_DTYPES_NP_TO_GDAL = {  np.dtype(np.float32):gdalconst.GDT_Float32,
                        np.dtype(np.float64):gdalconst.GDT_Float64,
                        np.dtype(np.int16):gdalconst.GDT_Int16,
                        np.dtype(np.int32):gdalconst.GDT_Int32,
                        np.dtype(np.uint16):gdalconst.GDT_UInt16,
                        np.dtype(np.int8):gdalconst.GDT_Byte,
                        np.dtype(np.uint8):gdalconst.GDT_UInt16}
_DTYPES_GDAL_TO_NP = {a_value:a_key for a_key,a_value in _DTYPES_NP_TO_GDAL.items()}

class OutsideExtent(Exception):
    pass

class RasterDataset(object):
    '''
    Encapsulates a basic 2-D GDAL raster object with
    common utility functions.
    '''

    def __init__(self, ds_path):
        '''
        Parameters
        ----------
        ds_path : str
            A file path to a raster dataset (eg GeoTIFF file)
        '''

        self.gdal_ds = gdal.Open(ds_path)
        self.ds_path = ds_path
        
        # GDAL GeoTransform.
        # Top left x,y are for the upper left corner of upper left pixel
        # GeoTransform[0] /* top left x */
        # GeoTransform[1] /* w-e pixel resolution */
        # GeoTransform[2] /* rotation, 0 if image is "north up" */
        # GeoTransform[3] /* top left y */
        # GeoTransform[4] /* rotation, 0 if image is "north up" */
        # GeoTransform[5] /* n-s pixel resolution */
        self.geo_t = np.array(self.gdal_ds.GetGeoTransform())

        self.projection = self.gdal_ds.GetProjection()
        self.source_sr = osr.SpatialReference()
        self.source_sr.ImportFromWkt(self.projection)
        self.target_sr = osr.SpatialReference()
        self.target_sr.ImportFromEPSG(_PROJECTION_GEO_WGS84)
        self.coordTrans_src_to_wgs84 = osr.CoordinateTransformation(self.source_sr,
                                                                    self.target_sr)
        self.coordTrans_wgs84_to_src = osr.CoordinateTransformation(self.target_sr,
                                                                    self.source_sr)

        self.min_x = self.geo_t[0]
        self.max_x = self.min_x + (self.gdal_ds.RasterXSize * self.geo_t[1])
        self.max_y = self.geo_t[3]
        self.min_y = self.max_y - (-self.gdal_ds.RasterYSize * self.geo_t[5])

        self.rows = self.gdal_ds.RasterYSize
        self.cols = self.gdal_ds.RasterXSize
        self.ndata = self.gdal_ds.GetRasterBand(1).GetNoDataValue()


    def get_coord_mesh_grid(self):
        '''
        Build a native projection coordinate mesh grid for the raster
        
        Returns
        ----------
        yGrid : ndarray
            A 2-D ndarray with the same shape as the raster containing
            the y-coordinate of the center of each grid cell.
        xGrid : ndarray
            A 2-D ndarray with the same shape as the raster containing
            the x-coordinate of the center of each grid cell.
        '''

        # Get the upper left and right point x coordinates in the raster's projection
        ulX = self.geo_t[0] + (self.geo_t[1] / 2.0)
        urX = self.get_coord(0, self.gdal_ds.RasterXSize - 1)[1]

        # Get the upper left and lower left y coordinates
        ulY = self.geo_t[3] + (self.geo_t[5] / 2.0)
        llY = self.get_coord(self.gdal_ds.RasterYSize - 1, 0)[0]

        # Build 1D arrays of x,y coords
        x = np.linspace(ulX, urX, self.gdal_ds.RasterXSize)
        y = np.linspace(ulY, llY, self.gdal_ds.RasterYSize)

        xGrid, yGrid = np.meshgrid(x, y)

        return yGrid, xGrid

    def get_coord_grid_1d(self):
        '''
        Build 1-D native projection coordinates for y and x dimensions
        of the raster.
        
        Returns
        ----------
        y : ndarray
            A 1-D array of the y-coordinates of each row.
        x : ndarray
            A 1-D array of the x-coordinates of each column.
        '''

        # Get the upper left and right point x coordinates in the raster's projection
        ulX = self.geo_t[0] + (self.geo_t[1] / 2.0)
        urX = self.get_coord(0, self.gdal_ds.RasterXSize - 1)[1]

        # Get the upper left and lower left y coordinates
        ulY = self.geo_t[3] + (self.geo_t[5] / 2.0)
        llY = self.get_coord(self.gdal_ds.RasterYSize - 1, 0)[0]

        # Build 1D arrays of x,y coords
        x = np.linspace(ulX, urX, self.gdal_ds.RasterXSize)
        y = np.linspace(ulY, llY, self.gdal_ds.RasterYSize)

        return y, x

    def get_coord(self, row, col):
        '''
        Get the native projection coordinates for a specific grid cell
        
        Parameters
        ----------
        row : int
            The row of the grid cell (zero-based)
        col : int
            The column of the grid cell (zero-based)
            
        Returns
        ----------
        yCoord : float
            The y projection coordinate of the grid cell.
        xCoord : float
            The x projection coordinate of the grid cell.   
        '''

        xCoord = ((self.geo_t[0] + col * self.geo_t[1] + row * self.geo_t[2]) + 
                  self.geo_t[1] / 2.0)
        yCoord = ((self.geo_t[3] + col * self.geo_t[4] + row * self.geo_t[5]) + 
                  self.geo_t[5] / 2.0)

        return yCoord, xCoord

    def get_row_col(self, lon, lat, check_bounds=True):
        '''
        Get the row, column grid cell offset for the raster based on an input
        WGS84 longitude, latitude point.
        
        Parameters
        ----------
        lon : float
            The longitude of the point
        lat : float
            The latitude of the point
        check_bounds : bool
            If True, will check if the point is within the bounds
            of the raster and raise a ValueError if not. If set to False and
            point is outside raster, the returned row, col will be clipped to
            the raster edge. 
        Returns
        ----------
        row : int
            The row of the closet grid cell to the lon,lat point (zero-based)
        col : int
            The column of the closet grid cell to the lon,lat point (zero-based) 
        '''

        xGeo, yGeo, zGeo = self.coordTrans_wgs84_to_src.TransformPoint(lon, lat)
        
        if check_bounds:

            if not self.is_inbounds(xGeo, yGeo):
                raise ValueError("lat/lon outside raster bounds: " + str(lat) + 
                                 "," + str(lon))

        originX = self.geo_t[0]
        originY = self.geo_t[3]
        pixelWidth = self.geo_t[1]
        pixelHeight = self.geo_t[5]

        xOffset = np.abs(np.int((xGeo - originX) / pixelWidth))
        yOffset = np.abs(np.int((yGeo - originY) / pixelHeight))
        
        # clip row,col if outside raster bounds
        row = self._check_cellxy_valid(int(yOffset), self.rows)
        col = self._check_cellxy_valid(int(xOffset), self.cols)
            
        return row, col

    def get_data_value(self, lon, lat):
        '''
        Get the nearest grid cell data value to an input
        WGS84 longitude, latitude point. Will raise an 
        OutsideExtent exception if the longitude, latitude 
        is outside the bounds of the raster.
        
        Parameters
        ----------
        lon : float
            The longitude of the point
        lat : float
            The latitude of the point
            
        Returns
        ----------
        data_val : dtype of raster
            The data value of the closet grid cell to the lon,lat point         
        '''

        row, col = self.get_row_col(lon, lat)
        data_val = self.gdal_ds.ReadAsArray(col, row, 1, 1)[0, 0]

        return data_val

    def read_as_array(self):
        '''
        Read in the entire raster as a 2-D array
        
        Returns
        ----------
        a : MaskedArray
            A 2-D MaskedArray of the raster data. No data values
            are masked.      
        '''

        a = self.gdal_ds.GetRasterBand(1).ReadAsArray()
        a = np.ma.masked_equal(a, self.gdal_ds.GetRasterBand(1).GetNoDataValue())
        return a
    
    def output_new_ds(self, fpath, a, gdal_dtype, ndata=None, gdal_driver="GTiff"):
        '''
        Output a new raster with same geotransform and projection as this RasterDataset
        
        Parameters
        ----------
        fpath : str
            The filepath for the new raster.
        a : ndarray or MaskedArray
            A 2-D array of the raster data to be output.
            If MaskedArray, masked values are set to ndata
        gdal_dtype : str
            A gdal datatype from gdalconst.GDT_*.
        ndata : num, optional
            The no data value for the raster
        gdal_driver : str, optional
            The GDAL driver for the output raster data format
        '''
        
        ds_out = gdal.GetDriverByName(gdal_driver).Create(fpath, int(a.shape[1]),
                                                          int(a.shape[0]), 1, gdal_dtype)
        ds_out.SetGeoTransform(self.gdal_ds.GetGeoTransform())
        ds_out.SetProjection(self.gdal_ds.GetProjection())
        
        band_out = ds_out.GetRasterBand(1)
        if ndata is not None:
            band_out.SetNoDataValue(ndata)
        band_out.WriteArray(np.ma.filled(a, ndata))
                    
        ds_out.FlushCache()
        ds_out = None

    def resample_to_ds(self, fpath, ds_src, gdal_gra, gdal_driver="GTiff"):
        '''
        Resample a different RasterDataset to this RasterDataset's grid
        
        Parameters
        ----------
        fpath : str
            The filepath for the new resampled raster.
        ds_src : RasterDataset
            The RasterDataset to  be resampled
        gdal_gra : str
            A gdal resampling algorithm from gdalconst.GRA_*.
        gdal_driver : str, optional
            The GDAL driver for the output raster data format
            
        Returns
        ----------
        grid_out : RasterDataset
            A RasterDataset pointing to the resampled raster  
        '''
        
        grid_src = ds_src.gdal_ds
        grid_dst = self.gdal_ds
        
        proj_src = grid_src.GetProjection()
        dtype_src = grid_src.GetRasterBand(1).DataType
        ndata_src = grid_src.GetRasterBand(1).GetNoDataValue()
        
        proj_dst = grid_dst.GetProjection()
        geot_dst = grid_dst.GetGeoTransform()
            
        grid_out = gdal.GetDriverByName(gdal_driver).Create(fpath, grid_dst.RasterXSize,
                                                            grid_dst.RasterYSize, 1, dtype_src)
        
        if ndata_src is not None:
            band = grid_out.GetRasterBand(1)
            band.Fill(ndata_src)
            band.SetNoDataValue(ndata_src)
        
        grid_out.SetGeoTransform(geot_dst)
        grid_out.SetProjection(proj_dst)
        grid_out.FlushCache()
        
        gdal.ReprojectImage(grid_src, grid_out, proj_src, proj_dst, gdal_gra)
        grid_out.FlushCache()
        # Make sure entire grid is written by setting to None. 
        # FlushCache doesn't seem to write the entire grid after resampling?
        grid_out = None
        
        # return as RasterDataset
        grid_out = RasterDataset(fpath)
        
        return grid_out
        
    def is_inbounds(self, x_geo, y_geo):
        return (x_geo >= self.min_x and x_geo <= self.max_x and
                y_geo >= self.min_y and y_geo <= self.max_y)
    
    def _check_cellxy_valid(self, i, n):
        
        if i < 0:
            i = 0
        elif i >= n:
            i = n - 1
        
        return i
        
def init_out_nc(fpath_out, lon, lat, dates, vname, vname_attrs, ds_attrs,
                xy_dims=None, time_bnds_deltas=None, **kwargs):
    ''' Initialize a netcdf file for a 3D variable: (time, y, x)
    
        Parameters
        ----------
        fpath_out : str
            The filepath for the new netcdf file
        lon : array
            Longitude values for x dimension
        lat : array
            Latitude values for y dimension
        dates : DatatimeIndex
            Datetimes for time dimension
        vname : str
            Name of variable
        vname_attrs : dict-like
            Dictionary of variable attributes
        ds_attrs : dict-like
            Dictionary of global dataset attributes
        xy_dims : tuple, optional
            A tuple of length 2 containing arrays of the x and y dimension values
            if the x and y dimensions are not just lon, lat
        time_bnds_deltas : tuple, optional
            A tuple of 2 pandas.Timedelta objects specifying the time deltas that
            should be subtracted and added to the time values to create
            a time_bnds variable
        **kwargs :
            Keywords arguments that are passed on to netCDF4.createVariable for
            creating the main 3D variable
            
        Returns
        ----------
        netCDF4.Dataset
            Dataset object pointing to new netcdf file
    '''
    
    ds_out = nc.Dataset(fpath_out, 'w')
    ds_out.setncatts(dict(ds_attrs))
    
    # Create time variable
    ds_out.createDimension('time', dates.size)
    times = ds_out.createVariable('time', 'f8', ('time',), fill_value=False)
    times.long_name = "time"
    times.units = "days since 1900-01-01 00:00:00"
    times.standard_name = "time"
    times.calendar = "standard"
    times[:] = nc.date2num(dates.to_pydatetime(), times.units)
    
    if time_bnds_deltas is not None:
        
        ds_out.createDimension('nv', 2)
        times.bounds = 'time_bnds'
        
        time_starts = nc.date2num((dates - time_bnds_deltas[0]).to_pydatetime(),
                                  units=times.units)
        time_ends = nc.date2num((dates + time_bnds_deltas[1]).to_pydatetime(),
                                units=times.units)
        a_time_bnds = np.array([time_starts,time_ends]).T
        
        time_bnds = ds_out.createVariable('time_bnds', 'f8', ('time', 'nv'),
                                          fill_value=False)
        time_bnds[:] = a_time_bnds
    
    
    if xy_dims is None:
    
        ds_out.createDimension('lon', lon.size)
        ds_out.createDimension('lat', lat.size)
        dims_lon = ('lon',)
        dims_lat = ('lat',)
        dims_main = ('time','lat','lon')
        
    else:
        
        ds_out.createDimension('x', xy_dims[0].size)
        ds_out.createDimension('y', xy_dims[1].size)
        
        # Create x variable
        x_var = ds_out.createVariable('x', 'f8', ('x',), fill_value=False)
        x_var.long_name = 'x coordinate of projection'
        x_var.standard_name = 'projection_x_coordinate'
        x_var.units = 'm'
        x_var[:] = xy_dims[0]
        
        # Create y variable
        y_var = ds_out.createVariable('y', 'f8', ('y',), fill_value=False)
        y_var.long_name = 'y coordinate of projection'
        y_var.standard_name = 'projection_y_coordinate'
        y_var.units = 'm'
        y_var[:] = xy_dims[1]
        
        dims_lon = ('y','x')
        dims_lat = ('y','x')
        dims_main = ('time','y','x')
        

    # Create longitude variable
    lon_var = ds_out.createVariable('lon', 'f8', dims_lon, fill_value=False)
    lon_var.long_name = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.units = 'degrees_east'
    lon_var[:] = lon
    
    # Create latitude variable
    lat_var = ds_out.createVariable('lat', 'f8', dims_lat, fill_value=False)
    lat_var.long_name = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.units = 'degrees_north'
    lat_var[:] = lat
    
    # Create main variable
    if not kwargs.has_key('datatype'):
        kwargs['datatype'] = 'f4'
    
    main_var = ds_out.createVariable(varname=vname, dimensions=dims_main, **kwargs)
    main_var.setncatts(dict(vname_attrs))
    ds_out.sync()
    
    return ds_out        

def open_as_da(fpath_raster):
    '''Open a georeferenced raster file as a DataArray
    
    Parameters
    ----------
    fpath_raster : str
        Filepath of the raster file
        
    Returns
    ----------
    xarray.DataArray
        Proj.4 projection string of raster file is specified by the "crs"
        attribute.
    '''
    
    
    ds = RasterDataset(fpath_raster)
        
    a = ds.read_as_array()
    
    if ds.source_sr.IsGeographic():
        
        lat,lon = ds.get_coord_grid_1d()
        da = xr.DataArray(a, dims=['lat','lon'], coords=[lat, lon])
        
    else: 
        
        y,x = ds.get_coord_grid_1d()
        y_grid, x_grid = ds.get_coord_mesh_grid()
        p1 = pyproj.Proj(ds.source_sr.ExportToProj4())
        p2 = pyproj.Proj('+init=epsg:4326 +no_defs')
        lon,lat = pyproj.transform(p1, p2, x_grid.ravel(), y_grid.ravel())
        lon = xr.DataArray(lon.reshape(x_grid.shape),dims=['y','x'],coords=[y,x])
        lat = xr.DataArray(lat.reshape(y_grid.shape),dims=['y','x'],coords=[y,x])
    
        da = xr.DataArray(a,dims=['y','x'],
                          coords={'y':y,'x':x,'lon':lon,'lat':lat})
    
    da.attrs['crs'] = ds.source_sr.ExportToProj4()
    
    return da
    

def to_raster_file(da, fpath_out, xname='lon', yname='lat', crs='+init=epsg:4326 +no_defs',
                   ndata=None, dtype=None, gdal_driver="GTiff"):
    '''Output a DataArray as a georeferenced raster file
    
    Parameters
    ----------
    da : DataArray
        A 2D DataArray
    fpath_out : str
        Filepath for the output raster file
    xname : str, optional
        Name of DataArray's x dimension. Default: lon
    yname : str, optional
        Name of DataArray's y dimension. Default: lat
    crs : str, optional
        proj.4 CRS string for the DataArray's projection. Default: WGS84
    ndata : numeric, optional
        The no data value for the raster. All NA values in the DataArray
        will be filled with this value before writing
    dtype : numpy.dtype. optional
        The data type of the raster. If not specified, will use datatype
        of the DataArray
    gdal_driver : str, optional
        The GDAL driver for the output raster data format. Default: GTiff
        
    '''

    dtype = np.dtype(da.dtype) if dtype is None else np.dtype(dtype)
    gdal_dtype = _DTYPES_NP_TO_GDAL[dtype]
    
    xvals = da[xname].values
    yvals = da[yname].values
    
    geo_t = [None]*6
    #n-s pixel height/resolution needs to be negative
    geo_t[5] = -np.abs(yvals[0] - yvals[1])   
    geo_t[1] = np.abs(xvals[0] - xvals[1])
    geo_t[2],geo_t[4] = (0.0,0.0)
    geo_t[0] = xvals[0] - (geo_t[1]/2.0) 
    geo_t[3] = yvals[0] + np.abs(geo_t[5]/2.0)
    
    ds_out = gdal.GetDriverByName(gdal_driver).Create(fpath_out, int(da.shape[1]),
                                                      int(da.shape[0]),
                                                      1, gdal_dtype)
    ds_out.SetGeoTransform(geo_t)
    
    sr = osr.SpatialReference()
    sr.ImportFromProj4(crs)
    ds_out.SetProjection(sr.ExportToWkt())
    
    band_out = ds_out.GetRasterBand(1)
    if ndata is not None:
        band_out.SetNoDataValue(ndata)
        da = da.fillna(ndata)
        
    band_out.WriteArray(da.values)
    ds_out.FlushCache()
    ds_out = None
    

def to_inmemory_rio(da, xname='lon', yname='lat', crs='+init=epsg:4326 +no_defs'):
    '''Create a rasterio MemoryFile from an xarray DataArray object
    
    Parameters
    ----------
    da : DataArray
        A 2D DataArray
    xname : str, optional
        Name of DataArray's x dimension. Default: lon
    yname : str, optional
        Name of DataArray's y dimension. Default: lat
    crs : str, optional
        proj.4 CRS string for the DataArray's projection. Default: WGS84
        
    Returns
    ----------
    MemoryFile
        The rasterio MemoryFile
    '''
    
    x_size = abs(da[xname].values[0] - da[xname].values[1])
    y_size = abs(da[yname].values[0] - da[yname].values[1])
    x_topleft = da[xname].values[0] - (x_size / 2.0)
    y_topleft = da[yname].values[0] + (y_size / 2.0)
    a_trans = rio.transform.from_origin(x_topleft, y_topleft,
                                        x_size, y_size)
    a_crs = rio.crs.CRS.from_string(crs)
    
    meta_out = {'count':1, 'crs': a_crs,
                'driver': 'GTiff', 'dtype': da.dtype,
                'height': da[yname].size, 'width': da[xname].size,
                'transform': a_trans, 'nodata':-9999}
    
    a_out = da.fillna(-9999).values
    a_out = a_out.reshape(1, *a_out.shape)
    
    with rio.MemoryFile() as memfile:
        ds_mem = memfile.open(**meta_out)
        ds_mem.write(a_out)
    
    return ds_mem


def mask_and_crop(da, features, xname='lon', yname='lat', crs='+init=epsg:4326 +no_defs'):
    '''Mask areas in DataArray outside shapes and crop to extent of shapes
    
    Parameters
    ----------
    da : DataArray
        A 2D DataArray
    features : list of polygons
        A list of polygon shapes
    xname : str, optional
        Name of DataArray's x dimension. Default: lon
    yname : str, optional
        Name of DataArray's y dimension. Default: lat
    crs : str, optional
        proj.4 CRS string for the DataArray's projection. Default: WGS84
     
    Returns
    ----------
    DataArray
        The cropped and masked DataArray
    
    '''

    ds_mem = to_inmemory_rio(da,xname, yname, crs)
    a_crop, a_transform = rio_mask.mask(ds_mem, features, nodata=-9999,
                                       crop=True, all_touched=True)
    a_crop = a_crop.data
    a_crop[a_crop == -9999] = np.nan
    
    out_meta = ds_mem.meta.copy()
    out_meta.update({"height": a_crop.shape[1],
                     "width": a_crop.shape[2],
                     "transform": a_transform})
    
    ds_mem.close()
    del ds_mem
    
    with rio.MemoryFile() as memfile:
        
        with memfile.open(**out_meta) as ds_mem:
        
            x_s = np.linspace(ds_mem.xy(0, 0)[0], ds_mem.xy(0, ds_mem.width - 1)[0],
                              num=ds_mem.width)
            y_s = np.linspace(ds_mem.xy(0, 0)[1], ds_mem.xy(ds_mem.height - 1, 0)[1],
                              num=ds_mem.height)
    
    if da.__contains__('time'):
    
        da_out = xr.DataArray(a_crop[0], coords={xname:x_s, yname:y_s,
                                                'time': da['time'].values},
                              dims=[yname, xname], name=da.name)
        da_out.time.attrs = da.time.attrs
        
    else:
        
        da_out = xr.DataArray(a_crop[0], coords={xname:x_s, yname:y_s},
                              dims=[yname, xname], name=da.name)
        
        
    da_out[yname].attrs = da[yname].attrs
    da_out[xname].attrs = da[xname].attrs
    da_out.attrs = da.attrs
    
    return da_out

def _get_transform(da, xname, yname):
    
    x_size = abs(da[xname].values[0] - da[xname].values[1])
    y_size = abs(da[yname].values[0] - da[yname].values[1])
    x_topleft = da[xname].values[0] - (x_size / 2.0)
    y_topleft = da[yname].values[0] + (y_size / 2.0)
    
    a_trans = rio.transform.from_origin(x_topleft, y_topleft, x_size, y_size)
    
    return a_trans
    
def reproject_da(da_src, da_dst, crs_src, crs_dst, xy_src, xy_dst, resampling='nearest'):
    '''Reproject a DataArray to a different DataArray grid
    
    Parameters
    ----------
    da_src : DataArray
        A 2D DataArray with the source data
    da_dst : DataArray
        A 2D DataArray representing the destination grid
    crs_src : str
        proj.4 CRS string for the source DataArray
    crs_dst : str
        proj.4 CRS string for the destination DataArray
    xy_src : tuple
        x, y dimension names for the source DataArray
    xy_dest : tuple
        x, y dimension names for the destination DataArray
    resampling : str
        The resampling method. One of: nearest, bilinear, cubic, cubic_spline,
        average, gauss, lanczos, max, med, min, mode, q1, q3. Default: nearest
     
    Returns
    ----------
    DataArray
        The reprojected DataArray
    '''
    
    src_transform = _get_transform(da_src,xy_src[0],xy_src[1])
    src_crs = rio.crs.CRS.from_string(crs_src)
    
    dst_transform = _get_transform(da_dst,xy_dst[0],xy_dst[1])
    dst_crs = rio.crs.CRS.from_string(crs_dst)
    da_dst = da_dst.copy()
    da_dst.attrs = da_src.attrs
    
    reproject(da_src.values,da_dst.values,src_transform=src_transform,
              src_crs=src_crs,src_nodata=np.nan, dst_transform=dst_transform,
              dst_crs=dst_crs,resampling=Resampling[resampling])
    
    return da_dst
    