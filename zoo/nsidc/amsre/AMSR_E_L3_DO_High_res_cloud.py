"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an NSIDC AMSR grid
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AMSR_E_L3_DO_High_res_cloud.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_GDAL = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'AMSR_E_L3_DailyOcean_V03_20020619.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'High_res_cloud'

    if USE_GDAL:

        import gdal
        import mpl_toolkits.basemap.pyproj as pyproj

        GRID_NAME = 'GlobalGrid'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)

        # Construct the grid.  The projection is GEO, so this immediately
        # gives us latitude and longitude.
        meta = gdset.GetMetadata()
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
        x = np.linspace(x0, x0 + xinc*nx, nx)
        y = np.linspace(y0, y0 + yinc*ny, ny)
        longitude, latitude = np.meshgrid(x, y)

        scale_factor = float(meta['Scale'])
        units = meta['Unit']
        del gdset

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.float64)

        # Retrieve attributes.
        attrs = data2D.attributes(full=1)
        scale_factor = attrs["Scale"][0]
        units = attrs["Unit"][0]

        # Read global attribute.
        gridmeta = hdf.attributes(full=1)["StructMetadata.0"][0]

        # Construct the grid.  The needed information is in a global attribute
        # called 'StructMetadata.0'.  Use regular expressions to tease out the
        # extents of the grid.  In addition, the grid is in packed decimal
        # degrees, so we need to normalize to degrees.

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
        x = np.linspace(x0, x1, nx)
        y = np.linspace(y0, y1, ny)
        longitude, latitude = np.meshgrid(x, y)

    # Apply the attributes information.
    data[data == -9999] = np.nan
    data = data * scale_factor
    data = np.ma.masked_array(data, np.isnan(data))

    long_name = DATAFIELD_NAME.replace('_', ' ')

    m = Basemap(projection='cyl', resolution='l', lon_0=0,
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 45), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
