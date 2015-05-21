"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LAADS MYD (MODIS-
AQUA) swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MYD07_L2_Water_Vapor.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = False


def run():

    # If a certain environment variable is set, look there for the input
    # files, otherwise look in the current directory.
    FILE_NAME = 'MYD07_L2.A2002184.2200.005.2006133121629.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'Water_Vapor'

    if USE_NETCDF4:

        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)

        # The netCDF4 module will correctly apply fill value and
        # scaling equation.
        data = nc.variables[DATAFIELD_NAME][:]

        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]

        # Retrieve attributes.
        units = nc.variables[DATAFIELD_NAME].units
        long_name = nc.variables[DATAFIELD_NAME].long_name

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.double)

        # Read geolocation dataset.
        lat = hdf.select('Latitude')
        latitude = lat[:]
        lon = hdf.select('Longitude')
        longitude = lon[:]

        # Retrieve attributes.
        attrs = data2D.attributes(full=1)
        long_name = attrs["long_name"][0]
        add_offset = attrs["add_offset"][0]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["scale_factor"][0]
        valid_min = attrs["valid_range"][0][0]
        valid_max = attrs["valid_range"][0][1]
        units = attrs["units"][0]

        # Apply _FillValue, scale and offset.
        invalid = np.logical_or(data > valid_max,
                                data < valid_min)
        invalid = np.logical_or(invalid, data == _FillValue)
        data[invalid] = np.nan
        data = data * scale_factor + add_offset
        data = np.ma.masked_array(data, np.isnan(data))

    # The data is local to Alaska, so no need for a global or hemispherical
    # projection.
    m = Basemap(projection='laea', resolution='l',
                lat_ts=65, lat_0=65, lon_0=-150,
                width=4800000, height=3500000)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(40, 81, 10), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-210, -89., 30), labels=[0, 0, 0, 1])
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
