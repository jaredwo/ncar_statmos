import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt

def map_pts(lon, lat, colors=None, cmap='viridis'):
    
    # plot in UTM Zone 17
    ax = plt.axes(projection=ccrs.UTM(17))
        
    shapename = 'admin_1_states_provinces_lakes_shp'
    states_shp = shpreader.natural_earth(resolution='50m',
                                         category='cultural', name=shapename)
    ax.add_geometries(shpreader.Reader(states_shp).geometries(),
                      ccrs.PlateCarree(),facecolor='none',lw=.5,edgecolor='k') 
    
    if colors is None:
        cmap = None
        
    sc = ax.scatter(lon, lat, c=colors, cmap=cmap, transform=ccrs.PlateCarree(),
                    edgecolor='k')
    
    ax.set_extent((lon.min(), lon.max(), lat.min(), lat.max()),crs=ccrs.Geodetic())

    if colors is not None:
        cb = plt.colorbar(sc, ax=ax)
    
    plt.show()