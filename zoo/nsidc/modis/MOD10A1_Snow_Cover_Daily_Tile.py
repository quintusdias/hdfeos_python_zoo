"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a NSIDC Level-2
MODIS Grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD10A1_Snow_Cover_Daily_Tile.py

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
    FILE_NAME = 'MOD10A1.A2000065.h00v08.005.2008237034422.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Snow_Cover_Daily_Tile'

    if USE_GDAL:    
        import gdal

        GRID_NAME = 'MOD_Grid_Snow_500m'
    
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)

        data = gdset.ReadAsArray()

        # Construct the grid.
        meta = gdset.GetMetadata()
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)


    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.float64)

        # Read global attribute.
        fattrs = hdf.attributes(full=1)
        gridmeta = fattrs["StructMetadata.0"][0]
            
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

    # In basemap, the sinusoidal projection is global, so we won't use it.
    # Instead we'll convert the grid back to lat/lons.
    sinu = pyproj.Proj("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
    wgs84 = pyproj.Proj("+init=EPSG:4326") 
    lon, lat= pyproj.transform(sinu, wgs84, xv, yv)

    # There's a wraparound issue for the longitude, as part of the tile extends
    # over the international dateline, and pyproj wraps longitude values west
    # of 180W (< -180) into positive territory.  Basemap's pcolormesh method
    # doesn't like that.
    lon[lon > 0] -= 360



    m = Basemap(projection='cyl', resolution='h',
                lon_0=-10,
                llcrnrlat=-5, urcrnrlat = 30,
                llcrnrlon=-185, urcrnrlon = -150)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(0, 21, 10), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, -159, 10), labels=[0, 0, 0, 1])

    # Use a discretized colormap since we have only four levels.
    # fill, ocean, no snow, missing
    # cmap = mpl.colors.ListedColormap(['black', 'blue', 'green', 'grey'])
    cmap = mpl.colors.ListedColormap(['grey',  'green', 'blue', 'black'])
    bounds = [0, 25, 39, 255, 256]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # 2400x2400 seems to be too much, so we'll subset it.
    m.pcolormesh(lon[::2,::2], lat[::2,::2], data[::2,::2], latlon=True,
                 cmap=cmap, norm=norm)
    
    color_bar = plt.colorbar()
    color_bar.set_ticks([12, 32, 147, 255.5])
    # color_bar.set_ticklabels(['fill', 'ocean', 'no snow', 'missing'])
    color_bar.set_ticklabels(['missing', 'no snow', 'ocean', 'fill'])
    color_bar.draw_all()

    basename = os.path.basename(FILE_NAME)
    long_name = 'Snow Cover Tile'
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
    
