"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a MOPITT HDF-EOS2
swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python TES_L2_O3_line_lvls.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

References
[1] http://tes.jpl.nasa.gov/uploadedfiles/TES_DPS_V11.8.pdf
"""

import datetime
import os

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'TES-Aura_L2-O3-Nadir_r0000011015_F05_07.he5'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    with h5py.File(FILE_NAME, mode='r') as f:

        group = '/HDFEOS/SWATHS/O3NadirSwath/Data Fields'
        data_var = f['/'.join([group, 'O3'])]
        o3_data = data_var[:]
        o3_longname = data_var.attrs['Title'].decode()
        o3_units = data_var.attrs['Units'].decode()
        o3_fillvalue = data_var.attrs['_FillValue']
        o3_missingvalue = data_var.attrs['MissingValue']
        o3_data[o3_data == o3_fillvalue] = np.nan
        o3_data[o3_data == o3_missingvalue] = np.nan

        pressure_var = f['/'.join([group, 'Pressure'])]
        pressure_data = pressure_var[:]
        pressure_longname = pressure_var.attrs['Title'].decode()
        pressure_units = pressure_var.attrs['Units'].decode()
        pressure_fillvalue = pressure_var.attrs['_FillValue']
        pressure_missingvalue = pressure_var.attrs['MissingValue']
        pressure_data[pressure_data == pressure_fillvalue] = np.nan
        pressure_data[pressure_data == pressure_missingvalue] = np.nan

        group = '/HDFEOS/SWATHS/O3NadirSwath/Geolocation Fields'
        time_var = f['/'.join([group, 'Time'])]
        time_data = time_var[:]

    # Time is second from TAI93.
    # See 4-25 of "TES Science Data Processing Standard and Special Observation
    # Data Products Specification" [1].
    # Please note that the computed time is off by 7 seconds from the
    # values stored in "/HDFEOS/SWATHS/O3NadirSwath/Data Fields/UTCTime".
    timebase = datetime.datetime(1993, 1, 1, 0, 0, 0)

    formatter = mpl.ticker.FormatStrFormatter('%.2g')
    basename = os.path.basename(FILE_NAME)

    xlabel = "{0} ({1})".format(o3_longname, o3_units)
    ylabel = "{0} ({1})".format(pressure_longname, pressure_units)

    ax1 = plt.subplot(2, 2, 1)
    ax1.semilogy(o3_data[55, :], pressure_data[55, :])
    ax1.set_ylabel(ylabel, fontsize=8)
    delta = datetime.timedelta(days=time_data[55] / 86400.0)
    timedatum = (timebase + delta).strftime('%d %a %Y %H:%M:%S')
    title = "{0}\n{1} at {2}".format(basename, o3_longname, timedatum)
    ax1.set_title(title, fontsize=8)
    ax1.set_xticks(np.arange(0e-6, 9e-6, 2e-6))
    ax1.xaxis.set_major_formatter(formatter)
    plt.tick_params(axis='both', labelsize=8)

    ax2 = plt.subplot(2, 2, 3)
    ax2.semilogy(o3_data[155, :], pressure_data[155, :])
    ax2.set_xlabel(xlabel, fontsize=8)
    ax2.set_ylabel(ylabel, fontsize=8)
    delta = datetime.timedelta(days=time_data[155] / 86400.0)
    timedatum = (timebase + delta).strftime('%d %a %Y %H:%M:%S')
    title = "{0}\n{1} at {2}".format(basename, o3_longname, timedatum)
    ax2.set_title(title, fontsize=8)
    ax2.xaxis.set_major_formatter(formatter)
    plt.tick_params(axis='both', labelsize=8)

    ax3 = plt.subplot(2, 2, 2)
    ax3.semilogy(o3_data[955, :], pressure_data[955, :])
    #    ax3.set_ylabel(ylabel, fontsize=8)
    delta = datetime.timedelta(days=time_data[955] / 86400.0)
    timedatum = (timebase + delta).strftime('%d %a %Y %H:%M:%S')
    title = "{0}\n{1} at {2}".format(basename, o3_longname, timedatum)
    ax3.set_title(title, fontsize=8)
    ax3.set_xticks(np.arange(0e-6, 9e-6, 2e-6))
    ax3.xaxis.set_major_formatter(formatter)
    plt.tick_params(axis='both', labelsize=8)

    ax4 = plt.subplot(2, 2, 4)
    ax4.semilogy(o3_data[1555, :], pressure_data[1555, :])
    ax4.set_xlabel(xlabel, fontsize=8)
    #   ax4.set_ylabel(ylabel, fontsize=8)
    delta = datetime.timedelta(days=time_data[1555] / 86400.0)
    timedatum = (timebase + delta).strftime('%d %a %Y %H:%M:%S')
    title = "{0}\n{1} at {2}".format(basename, o3_longname, timedatum)
    ax4.set_title(title, fontsize=8)
    ax4.xaxis.set_major_formatter(formatter)
    plt.tick_params(axis='both', labelsize=8)

    fig = plt.gcf()
    plt.show(block=False)

    pngfile = "{0}.py.l.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
