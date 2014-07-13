"""
This example code illustrates how to access and visualize a NSIDC
MODIS Grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD10C1_Day_CMG_Snow_Cover.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""

import os
import re

import gdal, osr
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

def run(FILE_NAME):
    
    # Identify the data field.
    GRID_NAME = 'MOD_CMG_Snow_5km'
    DATAFIELD_NAME = 'Day_CMG_Snow_Cover'
    
    gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                     GRID_NAME,
                                                     DATAFIELD_NAME)
    gdset = gdal.Open(gname)

    data = gdset.ReadAsArray()

    # Construct the grid.
    x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
    nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
    x = np.linspace(x0, x0 + xinc*nx, nx)
    y = np.linspace(y0, y0 + yinc*ny, ny)
    lon, lat = np.meshgrid(x, y)

    del gdset

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])

    # Bin the data as follows:
    # 0% snow
    # 1-99% snow
    # 100% snow
    # lake ice (107)
    # night (111)
    # cloud-obscured water (250)
    # water mask (254)
    # fill (255)
    lst = ['#00ff00', 
           '#888888',
           '#ffffff',
           '#ffafff',
           '#000000',
           '#63c6ff',
           '#0000cc',
           '#8928dd']
    cmap = mpl.colors.ListedColormap(lst)
    bounds = [0, 1, 100, 107, 111, 250, 254, 255, 256]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # Render the image in the projected coordinate system.
    m.pcolormesh(lon[::2,::2], lat[::2,::2], data[::2,::2],
                 latlon=True, cmap=cmap, norm=norm)
    
    color_bar = plt.colorbar()
    color_bar.set_ticks([0.5, 50, 103, 109, 180, 252, 254.5, 255.5])
    color_bar.set_ticklabels(['0% snow', '1-99% snow', '100% snow', 'lake ice',
                              'night', 'cloud-obscured water', 'water mask',
                              'fill'])
    color_bar.draw_all()
    fig = plt.gcf()
    
    plt.title('Day CMG Snow Cover')
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOD10C1.A2005018.005.2007349093349.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
