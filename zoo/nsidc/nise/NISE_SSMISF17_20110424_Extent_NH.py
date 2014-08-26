"""
This example code illustrates how to access and visualize a NSIDC NISE Grid file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python NISE_SSMISF17_20110424_Extent_NH.py

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

FILE_NAME = 'NISE_SSMISF17_20110424.HDFEOS'

def run(grid_file):
    
    # Identify the data field.
    GRID_NAME = 'Northern Hemisphere'
    DATAFIELD_NAME = 'Extent'
    
    gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(grid_file,
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

    # Reproject into WGS84
    lamaz = pyproj.Proj("+proj=laea +a=6371228 +lat_0=90 +lon_0=0 +units=m")
    wgs84 = pyproj.Proj("+init=EPSG:4326") 
    lon, lat= pyproj.transform(lamaz, wgs84, xv, yv)

    # Use a north polar azimuthal equal area projection.
    m = Basemap(projection='nplaea', resolution='l',
                boundinglat=40, lon_0=0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(0, 90, 15), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180, 45), labels=[0, 0, 0, 1])

    # Bin the data as follows:
    # 0 -- snow-free land
    # 1-20% sea ice -- blue
    # 21-40% sea ice -- blue-cyan
    # 41-60% sea ice -- blue
    # 61-80% sea ice -- cyan-blue
    # 81-100% sea ice -- cyan
    # 101 -- permanent ice
    # 103 -- dry snow
    # 252 mixed pixels at coastlines
    # 255 ocean
    lst = ['#004400', 
           '#0000ff',
           '#0044ff',
           '#0088ff',
           '#00ccff',
           '#00ffff',
           '#ffffff',
           '#440044',
           '#191919',
           '#000000',
           '#8888cc']
    cmap = mpl.colors.ListedColormap(lst)
    bounds = [0, 1, 21, 41, 61, 81, 101, 103, 104, 252, 255]
    tickpts = [0.5, 11, 31, 51, 71, 91, 102, 103.5, 178, 253.5] 
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # The corners cause trouble, so chop them out.
    idx = slice(5, 721)
    m.pcolormesh(lon[idx, idx], lat[idx, idx], data[idx, idx],
                 latlon=True, cmap=cmap, norm=norm)
    color_bar = plt.colorbar()
    color_bar.set_ticks(tickpts)
    color_bar.set_ticklabels(['snow-free\nland',
                              '1-20% sea ice',
                              '21-40% sea ice',
                              '41-60% sea ice',
                              '61-80% sea ice',
                              '81-100% sea ice',
                              'permanent\nice',
                              'dry\nsnow',
                              'mixed pixels\nat coastlines',
                              'ocean'])
    color_bar.draw_all()
    plt.title(DATAFIELD_NAME.replace('_',' '))

    fig = plt.gcf()
    plt.show()

    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.NH.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)

    del gdset


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    try:
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)
    except KeyError:
        pass

    run(FILE_NAME)
    
