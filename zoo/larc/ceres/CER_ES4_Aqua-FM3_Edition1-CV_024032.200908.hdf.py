"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LaRC CERES ES4 Aqua
 Grid HDF4 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python CER_ES4_Aqua-FM3_Edition1-CV_024032.200908.hdf.py

The HDF4 file must either be in your current working directory
or in a directory specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

# The file has many 'Net radiation flux' dataset under
# different Vgroup.
#
# Using netcdf4 will limit the number of datasets that you can plot.
#
# This is a good example why PyHDF is better than netcdf4 python.
USE_NETCDF4 = False

def run(FILE_NAME):

    # Identify the data field.
    DATAFIELD_NAME = 'Net radiant flux'

    if USE_NETCDF4:
        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)
    
        var = nc.variables[DATAFIELD_NAME]
        data = var[:].astype(np.float64)
        latitude = nc.variables['Colatitude'][:]
        longitude = nc.variables['Longitude'][:]

        # Read attributes.
        units = var.units
        fillvalue = var._FillValue
        long_name = var.long_name

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        # The file has many 'Net radiation flux' dataset under
        # different Vgroup.
        # 
        # Use HDFView to look up ref number.
        index = hdf.reftoindex(141)
        data1D = hdf.select(index)

        data = data1D[:].astype(np.double)

        # Read geolocation datasets.
        lat = hdf.select(hdf.reftoindex(185))
        latitude = lat[:]

        lon = hdf.select(hdf.reftoindex(184))
        longitude = lon[:]


        # Read attributes.
        attrs = data1D.attributes(full=1)
        ua=attrs["units"]
        units = ua[0]
        fva=attrs["_FillValue"]
        fillvalue = fva[0]
        lna=attrs["long_name"]
        long_name = lna[0]

    # Apply the fill value attribute.
    data[data == fillvalue] = np.nan
    data = np.ma.masked_array(data, np.isnan(data))

    # Adjust lat values.
    latitude = 90 - latitude

    plt.plot(latitude, data)
    plt.xlabel('Latitude (degrees_north)')
    plt.ylabel('{0} ({1})'.format(DATAFIELD_NAME, units))

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}\n{2} at Longitude={3} (degrees_east)'.format(basename, long_name, DATAFIELD_NAME, longitude[0]), fontsize=11)
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)
    
if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    ncfile = 'CER_ES4_Aqua-FM3_Edition1-CV_024032.200908.hdf'
    try:
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'], ncfile)
    except KeyError:
        fname = hdffile

    run(fname)
