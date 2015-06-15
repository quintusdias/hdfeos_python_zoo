"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC HIRDLS
HDF-EOS5 swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python HIRDLS_Aura_L2_v06_00_00_c01_2008d001.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import datetime
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = True


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'HIRDLS-Aura_L2_v06-00-00-c01_2008d001.he5'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        grp = nc.groups['HDFEOS'].groups['SWATHS'].groups['HIRDLS']
        data_var = grp.groups['Data Fields'].variables['O3']
        pres_var = grp.groups['Geolocation Fields'].variables['Pressure']
        time_var = grp.groups['Geolocation Fields'].variables['Time']

        # Read the data.
        data = data_var[0, :]
        pressure = pres_var[:]
        time = time_var[0]

        # Read the needed attributes.
        data_units = data_var.Units
        pres_units = pres_var.Units
        data_title = data_var.Title
        time_title = time_var.Title
        pres_title = pres_var.Title

    else:

        import h5py

        with h5py.File(FILE_NAME, mode='r') as f:

            dset_var = f['/HDFEOS/SWATHS/HIRDLS/Data Fields/O3']
            dset_pres = f['/HDFEOS/SWATHS/HIRDLS/Geolocation Fields/Pressure']
            dset_time = f['/HDFEOS/SWATHS/HIRDLS/Geolocation Fields/Time']

            # Read the data.
            data = dset_var[0, :]
            pressure = dset_pres[:]
            time = dset_time[0]

            # Read the needed attributes.
            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            data_units = dset_var.attrs['Units'].decode()
            pres_units = dset_pres.attrs['Units'].decode()
            data_title = dset_var.attrs['Title'].decode()
            time_title = dset_time.attrs['Title'].decode()
            pres_title = dset_pres.attrs['Title'].decode()

            fillvalue = dset_var.attrs['_FillValue']
            data[data == fillvalue] = np.nan

    # The date is stored as a six-digit number, YYYYMM.  Convert it into
    # a string.
    datestr = datetime.datetime(1993, 1, 1) + datetime.timedelta(seconds=time)

    # Apply log scale along the y-axis to get a better image.
    pressure = np.log10(pressure)
    plt.plot(data, pressure)

    # Save some screen space by using scientific notation for the xtick labels.
    formatter = plt.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-3, 4))
    plt.gca().xaxis.set_major_formatter(formatter)

    plt.xlabel('{0} ({1})'.format(data_title, data_units))
    plt.ylabel('{0} ({1})\nin log10 scale'.format(pres_title, pres_units))

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1} at {2}'.format(basename, data_title,
              datestr.strftime('%Y-%m-%d %H:%M:%S')))

    fig = plt.gcf()
    # plt.show()

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
