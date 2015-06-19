"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an LP DAAC MOD11C2
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD11C2_LST_Night_CMG.py

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
    FILE_NAME = 'MOD11C2.A2007073.005.2007098050130.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.

    DATAFIELD_NAME = 'LST_Night_CMG'

    if USE_GDAL:

        import gdal

        GRID_NAME = 'MODIS_8DAY_0.05DEG_CMG_LST'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)

        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)

        # Read the attributes.
        meta = gdset.GetMetadata()
        long_name = meta['long_name']
        units = meta['units']
        _FillValue = np.float(meta['_FillValue'])
        scale_factor = np.float(meta['scale_factor'])
        add_offset = np.float(meta['add_offset'])
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')]
        del gdset

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.double)

        # Read attributes.
        attrs = data2D.attributes(full=1)
        attrs = data2D.attributes(full=1)
        long_name = attrs["long_name"][0]
        valid_range = attrs["valid_range"][0]
        add_offset = attrs["add_offset"][0]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["scale_factor"][0]
        units = attrs["units"][0]

    # Have to be careful of the scaling equation here.
    invalid = data == _FillValue
    invalid = np.logical_or(invalid, data < valid_range[0])
    invalid = np.logical_or(invalid, data > valid_range[1])
    data[invalid] = np.nan
    data = (data - add_offset) * scale_factor
    data = np.ma.masked_array(data, np.isnan(data))

    # Normally we would use the HDF-EOS metadata to reconstruct the grid, but
    # the grid metadata is incorrect in this case, specifically the upper left
    # and lower right coordinates of the grid.  We'll construct the grid
    # manually.
    x = np.linspace(-180, 180, data.shape[1])
    y = np.linspace(90, -90, data.shape[0])
    lon, lat = np.meshgrid(x, y)

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 90, 45), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data)
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
