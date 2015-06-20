"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an LP DAAC MCD43C1
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MCD43C1_Black_Sky_Albedo_Num_Albedo_Bands1.py

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
    FILE_NAME = 'MCD43C1.A2006353.004.2007012185705.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Black_Sky_Albedo'

    if USE_GDAL:

        import gdal

        GRID_NAME = 'MOD_CMG_BRDF_0.05Deg'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)

        # Just retrieve the 2nd band.
        gdset = gdal.Open(gname)
        band = gdset.GetRasterBand(2)
        data = band.ReadAsArray().astype(np.float64)

        # Read attributes.
        meta = gdset.GetMetadata()
        long_name = meta['long_name']
        units = meta['units']
        _FillValue = np.float(meta['_FillValue'])
        scale_factor = np.float(meta['scale_factor'])
        add_offset = np.float(meta['add_offset'])
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')]
        dimname = "Num_Albedo_Bands"
        del gdset

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data3D = hdf.select(DATAFIELD_NAME)
        data = data3D[:, :, 1].astype(np.double)

        # Read attributes.
        attrs = data3D.attributes(full=1)
        long_name = attrs["long_name"][0]
        valid_range = attrs["valid_range"][0]
        add_offset = attrs["add_offset"][0]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["scale_factor"][0]
        units = attrs["units"][0]

        # Retrieve dimension name.
        dim = data3D.dim(2)
        dimname = dim.info()[0]

    invalid = data == _FillValue
    invalid = np.logical_or(invalid, data < valid_range[0])
    invalid = np.logical_or(invalid, data > valid_range[1])
    data[invalid] = np.nan
    data = scale_factor * (data - add_offset)
    data = np.ma.masked_array(data, np.isnan(data))

    # Normally we would use the following code to reconstruct the grid, but
    # the grid metadata is incorrect in this case, specifically the upper left
    # and lower right coordinates of the grid.  We'll construct the grid
    # manually, taking into account the fact that we're going to subset the
    # data by a factor of 10 (the grid size is 3600 x 7200).
    x = np.linspace(-180, 180, 720)
    y = np.linspace(90, -90, 360)
    lon, lat = np.meshgrid(x, y)

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 90, 45), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data[::10, ::10])
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}\nat {2}=1'.format(basename, long_name, dimname))

    fig = plt.gcf()
    plt.show(block=False)

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
