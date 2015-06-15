"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LaRC CERES AVG Grid
HDF4 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python CER_AVG_Aqua-FM3-MODIS_Edition2B_007005.200510.hdf.py

The netCDF file must either be in your current working directory
or in a directory specified by the environment variable HDFEOS_ZOO_DIR.
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
    FILE_NAME = 'CER_AVG_Aqua-FM3-MODIS_Edition2B_007005.200510.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'LW TOA Clear-Sky'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)

        # Subset the data to match the size of the swath geolocation fields.
        # Turn off autoscaling, we'll handle that ourselves due to presence of
        # a valid range.
        ncvar = nc.variables[DATAFIELD_NAME]
        data = ncvar[1, 0, :, :].astype(np.float64)
        latitude = nc.variables['Colatitude'][:]
        longitude = nc.variables['Longitude'][:]

        # Read attributes.
        units = ncvar.units
        fillvalue = ncvar._FillValue

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data4D = hdf.select(DATAFIELD_NAME)

        data = data4D[1, 0, :, :].astype(np.double)

        # Read geolocation datasets.
        lat = hdf.select('Colatitude')
        latitude = lat[:]
        lon = hdf.select('Longitude')
        longitude = lon[:]

        # Read attributes.
        attrs = data4D.attributes(full=1)
        units = attrs["units"][0]
        fillvalue = attrs["_FillValue"][0]

    # Apply the fill value attribute.
    data[data == fillvalue] = np.nan
    data = np.ma.masked_array(data, np.isnan(data))

    # Adjust lat/lon values.
    latitude = 90 - latitude
    longitude[longitude > 180] = longitude[longitude > 180] - 360

    # The data is global, so render in a global projection.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    long_name = DATAFIELD_NAME
    title = '{0}\n{1} at Monthly_Hourly_Avgs=0 and Stats=1'
    plt.title(title.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
