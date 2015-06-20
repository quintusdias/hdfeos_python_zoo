"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LAADS
MOD08 grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD08_D3_Cloud_Fraction_Liquid.py

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
    FILE_NAME = 'MOD08_D3.A2010001.005.2010006233008.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    GRID_NAME = 'mod08'
    DATAFIELD_NAME = 'Cloud_Fraction_Liquid'

    if USE_GDAL:

        import gdal

        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)

        # Read parameters for constructing the grid.
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)

        # Read fill value, valid range, scale factor, add_offset attributes.
        metadata = gdset.GetMetadata()
        valid_range = [float(x) for x in metadata['valid_range'].split(', ')]
        _FillValue = float(metadata['_FillValue'])
        scale_factor = float(metadata['scale_factor'])
        add_offset = float(metadata['add_offset'])
        units = metadata['units']
        long_name = metadata['long_name']

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.double)
        attrs = data2D.attributes(full=1)
        long_name = attrs["long_name"][0]
        add_offset = attrs["add_offset"][0]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["scale_factor"][0]
        valid_range = attrs["valid_range"][0]
        units = attrs["units"][0]

        # This product uses geographic projection.
        # Ideally, these parameters should be obtained by parsing
        # StructMetadta attribute but we assume that user has
        # checked it with HDFView.
        # Upper left corner: HDF-EOS2 convention
        x0 = -180
        y0 = 90
        # Grid spacing
        xinc = 1
        yinc = -1
        # Grid size
        nx = 360
        ny = 180

    invalid = data < valid_range[0]
    invalid = np.logical_or(invalid, data > valid_range[1])
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan

    data = scale_factor * (data - add_offset)
    data = np.ma.masked_array(data, np.isnan(data))

    x = np.linspace(x0, x0 + xinc*nx, nx)
    y = np.linspace(y0, y0 + yinc*ny, ny)
    lon, lat = np.meshgrid(x, y)

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}\n'.format(basename, long_name), fontsize=11)

    fig = plt.gcf()
    plt.show(block=False)

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
