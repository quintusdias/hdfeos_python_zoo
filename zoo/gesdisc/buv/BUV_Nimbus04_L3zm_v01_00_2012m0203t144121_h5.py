"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC MEaSUREs
Ozone Zonal Average HDF5 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python BUV_Nimbus04_L3zm_v01_00_2012m0203t144121_h5.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import datetime
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'BUV-Nimbus04_L3zm_v01-00-2012m0203t144121.h5'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        grp = nc.groups['Data_Fields']
        data_var = grp.variables['ProfileOzone']
        lat_var = grp.variables['Latitude']
        lev_var = grp.variables['ProfilePressureLevels']
        date_var = grp.variables['Date']

        # Read the data.
        data = data_var[0, :, :]
        lat = lat_var[:]
        lev = lev_var[:]
        date = date_var[0]

        # Read the needed attributes.
        data_units = data_var.units
        lat_units = lat_var.units
        lev_units = lev_var.units
        data_longname = data_var.long_name
        lat_longname = lat_var.long_name
        lev_longname = lev_var.long_name

    else:

        import h5py

        with h5py.File(FILE_NAME, mode='r') as f:

            dset_var = f['/Data_Fields/ProfileOzone']
            dset_lat = f['/Data_Fields/Latitude']
            dset_lev = f['/Data_Fields/ProfilePressureLevels']
            dset_date = f['/Data_Fields/Date']

            # Read the data.
            data = dset_var[0, :, :]
            lat = dset_lat[:]
            lev = dset_lev[:]
            date = dset_date[0]

            # Read the needed attributes.
            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            data_units = dset_var.attrs['units'].decode()
            lat_units = dset_lat.attrs['units'].decode()
            lev_units = dset_lev.attrs['units'].decode()
            data_longname = dset_var.attrs['long_name'].decode()
            lat_longname = dset_lat.attrs['long_name'].decode()
            lev_longname = dset_lev.attrs['long_name'].decode()

            # H5PY doesn't automatically turn the data into a masked array.
            fillvalue = dset_var.attrs['_FillValue']
            data[data == fillvalue] = np.nan
            data = np.ma.masked_array(data, np.isnan(data))

    # The date is stored as a six-digit number, YYYYMM.  Convert it into
    # a string.
    datestr = datetime.date(int(str(date)[0:4]), int(str(date)[4:6]), 1)

    # Apply log scale along the y-axis to get a better image.
    lev = np.log10(lev)
    plt.contourf(lat, lev, data.T)
    cb = plt.colorbar()
    cb.set_label('Unit:'+data_units)
    plt.xlabel('{0} ({1})'.format(lat_longname, lat_units))
    plt.ylabel('{0} ({1})\nin log10 scale'.format(lev_longname, lev_units))
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, data_longname))

    # Text position (80.0, -1.2) is relative to axes values.
    plt.text(80.0, -1.2, 'Date:{0}'.format(datestr.strftime('%Y-%m')),
             fontsize=12)
    fig = plt.gcf()
    plt.show(block=False)

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
