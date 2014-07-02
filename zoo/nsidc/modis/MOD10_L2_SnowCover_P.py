"""
This example code illustrates how to access and visualize a NSIDC Level-2
MODIS Swath data file in Python.

If you have any questions, suggestions, comments  on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD10_L2_SnowCover_P.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):
    
    dset = Dataset(FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Snow_Cover'
    
    # Subset the data to match the size of the swath geolocation fields.
    rows = slice(5, 4060, 10)
    cols = slice(5, 2708, 10)
    data = dset.variables['Snow_Cover'][rows, cols].astype(np.float64)
    latitude = dset.variables['Latitude'][:]
    longitude = dset.variables['Longitude'][:]
    
    # Draw a polar stereographic projection using the low resolution coastline
    # database.
    m = Basemap(projection='npstere', resolution='l',
                boundinglat=64, lon_0 = 0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(60.,81,10.))
    m.drawmeridians(np.arange(-180.,181.,30.), labels=[True,False,False,True])
    
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
    
    plt.title('Snow Cover')
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = basename + ".png"
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOD10_L2.A2000065.0040.005.2008235221207.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
