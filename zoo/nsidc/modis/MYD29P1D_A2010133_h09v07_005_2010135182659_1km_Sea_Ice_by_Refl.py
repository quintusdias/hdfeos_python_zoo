"""
This example code illustrates how to access and visualize an NSIDC MYD29
MODIS-AQUA 1km LAMAZ grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MYD29P1D_A2010133_h09v07_005_2010135182659_1km_Sea_Ice_by_Refl.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os
import re

import gdal
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np
from pyhdfeos.gd import GridFile
from pyhdf.SD import *


USE_GDAL = False

def run(FILE_NAME):
    
    GRID_NAME = 'MOD_Grid_Seaice_1km'
    DATAFIELD_NAME = 'Sea_Ice_by_Reflectance'
    
    if USE_GDAL:

        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray()

        meta = gdset.GetMetadata()
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
        x = np.linspace(x0, x0 + xinc*nx, nx)
        y = np.linspace(y0, y0 + yinc*ny, ny)
        xv, yv = np.meshgrid(x, y)

        # Reproject the coordinates out of lamaz into lat/lon.
        lamaz = pyproj.Proj("+proj=laea +a=6371228 +lat_0=90 +lon_0=0 +units=m")
        wgs84 = pyproj.Proj("+init=EPSG:4326") 
        lon, lat= pyproj.transform(lamaz, wgs84, xv, yv)

        del gdset

    else:

        gdf = GridFile(FILE_NAME)
        lat, lon = gdf.grids[GRID_NAME][:]

        sd_id = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        sds_id = sd_id.select(DATAFIELD_NAME)
        data = sds_id[:]

    # Draw a lambert equal area azimuthal basemap.
    m = Basemap(projection='laea', resolution='l', lat_ts=70,
                lat_0=70, lon_0=-180,
                width=2500000,height=2500000)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(50, 91, 10), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-220, -139, 10), labels=[0, 0, 0, 1])

    # Use a discretized colormap since we have only a few levels.
    # 0=missing data
    # 1=no decision
    # 11=night
    # 25=land
    # 37=inland water
    # 39=ocean
    # 50=cloud
    # 200=sea ice
    # 253=no input tile expected    
    # 254=non-production mask"
    # 255=fill
    lst = ['#727272',
           '#b7b7b7',
           '#ffff96',
           '#00ff00',
           '#232375',
           '#232375',
           '#63c6ff',
           '#ff0000',
           '#3f3f3f',
           '#000000',
           '#000000']
    cmap = mpl.colors.ListedColormap(lst)
    bounds = [0, 1, 11, 25, 37, 39, 50, 200, 253, 254, 255, 256]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    m.pcolormesh(lon, lat, data, latlon=True, cmap=cmap, norm=norm)
    color_bar = plt.colorbar()
    color_bar.set_ticks([0.5, 5.5, 18, 31, 38, 44.5, 125, 226.5, 253.5, 254.5, 255.5])
    color_bar.set_ticklabels(['missing', 'no decision', 'night', 'land',
                              'inland water', 'ocean', 'cloud', 'sea ice',
                              'no input tile expected', 'non-production mask',
                              'fill'])
    color_bar.draw_all()
    plt.title(DATAFIELD_NAME.replace('_', ' '))

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MYD29P1D.A2010133.h09v07.005.2010135182659.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
