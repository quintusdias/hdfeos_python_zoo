"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a MEaSURES VIP
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python VIP01P4_A2010001_002.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

In order for the netCDF code path to work, the netcdf library must be compiled
with HDF4 support.  Please see the README for details.
"""

import os
import re


import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

USE_NETCDF4 = False
USE_GDAL = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'VIP01P4.A2010001.002.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'CMG 0.05 Deg NDVI'

    if USE_GDAL:

        import gdal

        GRID_NAME = 'VIP_CMG_GRID'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)

        # Scale down the data by a factor of 6 so that low-memory machines
        # can handle it.
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)[::6, ::6]

        # Get any needed attributes.
        meta = gdset.GetMetadata()
        scale = np.float(meta['scale_factor'])
        fillvalue = np.float(meta['_FillValue'])
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')]
        units = meta['units']
        long_name = meta['long_name']

        # Construct the grid.
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        ny, nx = (gdset.RasterYSize / 6, gdset.RasterXSize / 6)
        x = np.linspace(x0, x0 + xinc*6*nx, nx)
        y = np.linspace(y0, y0 + yinc*6*ny, ny)
        lon, lat = np.meshgrid(x, y)

        del gdset

    else:

        if USE_NETCDF4:

            from netCDF4 import Dataset

            # The scaling equation isn't what netcdf4 expects, so turn it off.
            # Scale down the data by a factor of 6 so that low-memory machines
            # can handle it.
            nc = Dataset(FILE_NAME)
            ncvar = nc.variables[DATAFIELD_NAME]
            ncvar.set_auto_maskandscale(False)

            # Scale down the data by a factor of 6 so that low-memory machines
            # can handle it.
            data = ncvar[::6, ::6].astype(np.float64)

            # Get any needed attributes.
            scale = ncvar.scale_factor
            fillvalue = ncvar._FillValue
            valid_range = [np.float64(x)
                           for x in ncvar.valid_range.split(', ')]
            units = ncvar.units
            long_name = ncvar.long_name
            gridmeta = getattr(nc, 'StructMetadata.0')

        else:

            from pyhdf.SD import SD, SDC

            hdf = SD(FILE_NAME, SDC.READ)

            # Read dataset.
            data2D = hdf.select(DATAFIELD_NAME)
            data = data2D[:].astype(np.double)

            # Scale down the data by a factor of 6 so that low-memory machines
            # can handle it.
            data = data[::6, ::6]

            # Read attributes.
            attrs = data2D.attributes(full=1)
            long_name = attrs["long_name"][0]
            valid_range = attrs["valid_range"][0]
            valid_range = [np.float64(x) for x in valid_range.split(', ')]
            fillvalue = attrs["_FillValue"][0]
            scale = attrs["scale_factor"][0]
            units = attrs["units"][0]
            fattrs = hdf.attributes(full=1)
            gridmeta = fattrs["StructMetadata.0"][0]

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
        lon, lat = np.meshgrid(x, y)

    # Apply the attributes to the data.
    invalid = np.logical_or(data < valid_range[0], data > valid_range[1])
    invalid = np.logical_or(invalid, data == fillvalue)
    data[invalid] = np.nan
    data = data / scale
    data = np.ma.masked_array(data, np.isnan(data))

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])

    m.pcolormesh(lon, lat, data, latlon=True)
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
