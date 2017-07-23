import obsio
import pandas as pd

if __name__ == '__main__':
    
    path_local_ghcnd = '/storage/home/jwo118/group_store/topowx-2016/station_data'
    fpath_out = '/storage/home/jwo118/scratch/ncar_statmos/data/prcp_ghcnd_colorado_19480101_20151231.nc'
    download_ghcnd = False
    obs_elems = ('prcp', 'tobs_prcp')
    bbox = obsio.BBox(east_lon=-103.0, west_lon=-110.0, south_lat=36.0, north_lat=42.0)
    start_date = pd.Timestamp('1948-01-01')
    end_date = pd.Timestamp('2015-12-31')
    
    obsiof = obsio.ObsIoFactory(obs_elems, bbox, start_date, end_date)
    ghcnd = obsiof.create_obsio_dly_ghcnd(nprocs=3, bulk=True,
                                          local_data_path=path_local_ghcnd,
                                          download_updates=download_ghcnd)
    stns = ghcnd.stns
    
    ghcnd.to_netcdf(fpath_out, stns.station_id, start_date, end_date,
                    chk_rw=100, verbose=True)