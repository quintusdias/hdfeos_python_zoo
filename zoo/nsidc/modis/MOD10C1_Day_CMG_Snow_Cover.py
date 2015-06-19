"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

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
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

USE_GDAL = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'MOD10C1.A2005018.005.2007349093349.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Day_CMG_Snow_Cover'

    if USE_GDAL:

        import gdal

        GRID_NAME = 'MOD_CMG_Snow_5km'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray()

        # Read projection parameters.
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)

        del gdset

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.float64)

        # Read global attribute.
        fattrs = hdf.attributes(full=1)
        ga = fattrs["StructMetadata.0"]
        gridmeta = ga[0]

        # Read projection parameters.
        # The needed information is in a global attribute called
        # 'StructMetadata.0'.  Use regular expressions to tease out the extents
        # of the grid.
        ul_regex = re.compile(r'''UpperLeftPointMtrs=\(
                                  (?P<upper_left_x>[+-]?\d+\.\d+)
                                  ,
                                  (?P<upper_left_y>[+-]?\d+\.\d+)
                                  \)''', re.VERBOSE)
        match = ul_regex.search(gridmeta)
        x0 = np.float(match.group('upper_left_x')) / 1e6
        y0 = np.float(match.group('upper_left_y')) / 1e6

        lr_regex = re.compile(r'''LowerRightMtrs=\(
                                  (?P<lower_right_x>[+-]?\d+\.\d+)
                                  ,
                                  (?P<lower_right_y>[+-]?\d+\.\d+)
                                  \)''', re.VERBOSE)
        match = lr_regex.search(gridmeta)
        x1 = np.float(match.group('lower_right_x')) / 1e6
        y1 = np.float(match.group('lower_right_y')) / 1e6
        ny, nx = data.shape
        xinc = (x1 - x0) / nx
        yinc = (y1 - y0) / ny

    # Construct the grid.  It's already in lat/lon.
    x = np.linspace(x0, x0 + xinc*nx, nx)
    y = np.linspace(y0, y0 + yinc*ny, ny)
    lon, lat = np.meshgrid(x, y)

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
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
    m.pcolormesh(lon[::2, ::2], lat[::2, ::2], data[::2, ::2],
                 latlon=True, cmap=cmap, norm=norm)

    long_name = 'Day CMG Snow Cover'
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()

    color_bar = plt.colorbar(orientation='horizontal')
    color_bar.set_ticks([0.5, 50, 103, 109, 180, 252, 254.5, 255.5])
    color_bar.set_ticklabels(['0%\nsnow', '1-99%\nsnow', '100%\nsnow',
                              'lake\nice', 'night', 'cloud\n-obscured\nwater',
                              'water\nmask', 'fill'])
    color_bar.draw_all()

    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
