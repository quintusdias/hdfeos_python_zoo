"""
This example code illustrates how to access and visualize GESDISC_TRMM file in
Python.  If you have any questions, suggestions, comments  on this example,
please use the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to
see an  example of any other NASA HDF/HDF-EOS data product that is not 
listed in the HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo),
feel free to contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):

    DATAFIELD_NAME = 'surfaceRain'
    
    # Retrieve the data.
    dset = Dataset(FILE_NAME)
    data = dset.variables[DATAFIELD_NAME][:].astype(np.float64)
    units = dset.variables[DATAFIELD_NAME].units
    
    # Construct an indexed version of the data.
    levels = [0.0, 0.1, 1.0, 10.0, 30.0]
    Z = np.zeros(data.shape, dtype=np.float64)
    for j in range(len(levels)-1):
        Z[np.logical_and(data >= levels[j], data < levels[j+1])] = j  
    Z[data >= levels[-1]] = len(levels)
    
    
    # Retrieve the geolocation data.
    latitude = dset.variables['Latitude'][:]
    longitude = dset.variables['Longitude'][:]
    
    # There is a wrap-around effect to deal with.  Adjust the longitude by
    # modulus 360 to avoid the swath being smeared.
    longitude[longitude < -165] += 360
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-165, urcrnrlon = 197)
    
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 90., 30.))
    m.drawmeridians(np.arange(-180, 180., 45.))
    
    # Render the image in the projected coordinate system.
    x, y = m(longitude, latitude)
    
    # Use a discretized colormap since we have only five levels.
    colors = ['#0000ff', '#0088ff', '#8888ff', '#ff8888', '#ff0000']
    cmap = mpl.colors.ListedColormap(colors)
    bounds = np.arange(6)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    m.pcolormesh(x, y, Z, cmap=cmap, norm=norm)
    color_bar = plt.colorbar()
    color_bar.set_ticks([0.5, 1.5, 2.5, 3.5, 4.5])
    color_bar.set_ticklabels(['0', '0.1', '1.0', '10', '30'])
    
    fig = plt.gcf()
    
    plt.title('{0} ({1})'.format(DATAFIELD_NAME, units))
    plt.show()
    
    pngfile = "{0}.{1}.png".format(os.path.basename(FILE_NAME[:-4]),
                                   DATAFIELD_NAME)
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = '2A12.20140308.92894.7.HDF'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
