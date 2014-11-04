"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a NSIDC NISE Grid file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python NISE_SSMISF17_20110424_Extent_SH.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

USE_GDAL = False


def run(FILE_NAME):
    
    # Identify the data field.
    DATAFIELD_NAME = 'Extent'

    if USE_GDAL:
        import gdal
        GRID_NAME = 'Southern Hemisphere'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray()

        meta = gdset.GetMetadata()
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
        del gdset

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset. Dataset name 'Extent' exists under different groups.
        # Use reference number to resolve ambiguity.
        data2D = hdf.select(hdf.reftoindex(12))
        data = data2D[:,:].astype(np.float64)

        # Read global attribute.
        fattrs = hdf.attributes(full=1)
        ga = fattrs["StructMetadata.0"]
        gridmeta = ga[0]

        # Construct the grid.  The needed information is in a global attribute
        # called 'StructMetadata.0'.  Use regular expressions to tease out the
        # extents of the grid. 
        ul_regex = re.compile(r'''UpperLeftPointMtrs=\(
                                  (?P<upper_left_x>[+-]?\d+\.\d+)
                                  ,
                                  (?P<upper_left_y>[+-]?\d+\.\d+)
                                  \)''', re.VERBOSE)
        match = ul_regex.search(gridmeta)
        x0 = np.float(match.group('upper_left_x')) 
        y0 = np.float(match.group('upper_left_y')) 

        lr_regex = re.compile(r'''LowerRightMtrs=\(
                                  (?P<lower_right_x>[+-]?\d+\.\d+)
                                  ,
                                  (?P<lower_right_y>[+-]?\d+\.\d+)
                                  \)''', re.VERBOSE)
        match = lr_regex.search(gridmeta)
        x1 = np.float(match.group('lower_right_x'))
        y1 = np.float(match.group('lower_right_y'))
        ny, nx = data.shape
        xinc = (x1 - x0) / nx
        yinc = (y1 - y0) / ny

    x = np.linspace(x0, x0 + xinc*nx, nx)
    y = np.linspace(y0, y0 + yinc*ny, ny)
    xv, yv = np.meshgrid(x, y)

    # Reproject into WGS84
    lamaz = pyproj.Proj("+proj=laea +a=6371228 +lat_0=-90 +lon_0=0 +units=m")
    wgs84 = pyproj.Proj("+init=EPSG:4326") 
    lon, lat= pyproj.transform(lamaz, wgs84, xv, yv)

    # Use a south polar azimuthal equal area projection.
    m = Basemap(projection='splaea', resolution='l',
                boundinglat=-60, lon_0=0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 0, 15), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180, 30), labels=[0, 0, 0, 1])

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

    basename = os.path.basename(FILE_NAME)
    long_name = DATAFIELD_NAME
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.1.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'NISE_SSMISF17_20110424.HDFEOS'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
