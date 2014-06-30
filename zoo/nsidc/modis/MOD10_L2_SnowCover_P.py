"""
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):
    # The swath name is 'MOD_Swath_Snow', but we don't need that.
    
    # Identify the data field.
    DATAFIELD_NAME = 'Snow_Cover'
    
    rows = slice(5, 4060, 10)
    cols = slice(5, 2708, 10)
    
    dset = Dataset(FILE_NAME)
    data = dset.variables['Snow_Cover'][5:4060:10, 5:2708:10].astype(np.float64)
    
    latitude = dset.variables['Latitude'][:]
    longitude = dset.variables['Longitude'][:]
    
    # Draw a polar stereographic projection using the low resolution coastline
    # database.
    m = Basemap(projection='npstere', resolution='l',
                boundinglat=64, lon_0 = 0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(60.,81,10.))
    m.drawmeridians(np.arange(-180.,181.,30.), labels=[True,False,True,True])
    
    # Use a discretized colormap since we have only two levels.
    cmap = mpl.colors.ListedColormap(['grey','mediumblue'])
    bounds = [0, 19.5, 39]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # Render the image in the projected coordinate system.
    x, y = m(longitude, latitude)
    m.pcolor(x, y, data, alpha=0.90, cmap=cmap, norm=norm)
    
    # Must reset the alpha level to opaque for the colorbar.
    # See http://stackoverflow.com/questions/4478725/...
    # .../partially-transparent-scatter-plot-but-with-a-solid-color-bar
    color_bar = plt.colorbar()
    color_bar.set_alpha(1)
    color_bar.set_ticks([9.75, 29.25])
    color_bar.set_ticklabels(['missing data', 'ocean'])
    color_bar.draw_all()
    fig = plt.gcf()
    
    plt.show()
    plt.title('{0}\nSnow Cover'.format(os.path.basename(FILE_NAME)))
    
    fig.savefig('MOD10_L2.A2000065.0040.005.2008235221207_Snow_Cover_P_py.png')


if __name__ == "__main__":

    try:
        # If a certain environment variable is set, look there for the input
        # file.
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                             'MOD10_L2.A2000065.0040.005.2008235221207.hdf')
    except KeyError:
        # Look in the working directory.
        fname = 'MOD10_L2.A2000065.0040.005.2008235221207.hdf'

    run(fname)
    
