"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a NSIDC MOD29
Level 2 HDF-EOS2 Swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD29_A2013196_1250_005_2013196195940_hdf.py

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
    # file, otherwise look in the current directory.
    FILE_NAME = 'MOD29.A2013196.1250.005.2013196195940.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Ice_Surface_Temperature'

    if USE_NETCDF4:
        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)

        # Subset the data to match the size of the swath geolocation fields.
        # Turn off autoscaling, we'll handle that ourselves due to presence of
        # a valid range.
        var = nc.variables[DATAFIELD_NAME]
        var.set_auto_maskandscale(False)
        rows = slice(2, 2030, 5)
        cols = slice(2, 1354, 5)
        data = var[rows, cols].astype(np.float64)
        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]
        scale_factor = var.scale_factor
        add_offset = var.add_offset
        _FillValue = var._FillValue
        units = var.units
        valid_range = var.valid_range

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.float64)

        # Retrieve attributes.
        attrs = data2D.attributes(full=1)
        add_offset = attrs["add_offset"][0]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["scale_factor"][0]
        units = attrs["units"][0]
        valid_range = attrs["valid_range"][0]

        # Read lat and lon data from the matching geo-location file.
        GEO_FILE_NAME = 'MOD03.A2013196.1250.005.2013196194144.hdf'
        if 'HDFEOS_ZOO_DIR' in os.environ.keys():
            GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                         GEO_FILE_NAME)
        hdf_geo = SD(GEO_FILE_NAME, SDC.READ)

        # Read geolocation dataset from MOD03 product.
        lat = hdf_geo.select('Latitude')
        latitude = lat[:]
        lon = hdf_geo.select('Longitude')
        longitude = lon[:]

    # Apply the attributes.
    invalid = np.logical_or(data < valid_range[0],
                            data > valid_range[1])
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = data * scale_factor + add_offset
    datam = np.ma.masked_array(data, mask=np.isnan(data))

    # Draw a southern polar stereographic projection using the low resolution
    # coastline database.
    m = Basemap(projection='spstere', resolution='l',
                boundinglat=-64, lon_0=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-80, -59, 10))
    m.drawmeridians(np.arange(-180, 179, 30),
                    labels=[True, False, False, True])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    long_name = DATAFIELD_NAME
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    plt.show(block=False)
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
