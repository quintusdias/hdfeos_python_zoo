"""
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

def run(FILE_NAME):
    
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

    #fig = plt.figure(figsize = (15, 6))
    xlabel = "{0} ({1})".format(o3_longname, o3_units)
    ylabel = "{0} ({1})".format(pressure_longname, pressure_units)

    ax1 = plt.subplot(2, 2, 1)
    ax1.semilogy(o3_data[55,:], pressure_data[55,:])
    #ax1.set_xticks(
    #ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    delta =  datetime.timedelta(days=time_data[55]/86400.0)
    timedatum = (timebase + delta).strftime('%d %a %Y %H:%M:%S')
    ax1.set_title("{0} at {1}".format(o3_longname, timedatum), fontsize=12)
    ax1.set_xticks(np.arange(0e-6, 9e-6, 2e-6))
    ax1.xaxis.set_major_formatter(formatter)

    ax2 = plt.subplot(2, 2, 3)
    ax2.semilogy(o3_data[155,:], pressure_data[155,:])
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel)
    delta =  datetime.timedelta(days=time_data[155]/86400.0)
    timedatum = (timebase + delta).strftime('%d %a %Y %H:%M:%S')
    ax2.set_title("{0} at {1}".format(o3_longname, timedatum), fontsize=12)
    ax2.xaxis.set_major_formatter(formatter)

    ax3 = plt.subplot(2, 2, 2)
    ax3.semilogy(o3_data[955,:], pressure_data[955,:])
    ax3.set_ylabel(ylabel)
    delta =  datetime.timedelta(days=time_data[955]/86400.0)
    timedatum = (timebase + delta).strftime('%d %a %Y %H:%M:%S')
    ax3.set_title("{0} at {1}".format(o3_longname, timedatum), fontsize=12)
    ax3.set_xticks(np.arange(0e-6, 9e-6, 2e-6))
    ax3.xaxis.set_major_formatter(formatter)

    ax4 = plt.subplot(2, 2, 4)
    ax4.semilogy(o3_data[1555,:], pressure_data[1555,:])
    ax4.set_xlabel(xlabel)
    ax4.set_ylabel(ylabel)
    delta =  datetime.timedelta(days=time_data[1555]/86400.0)
    timedatum = (timebase + delta).strftime('%d %a %Y %H:%M:%S')
    ax4.set_title("{0} at {1}".format(o3_longname, timedatum), fontsize=12)
    ax4.xaxis.set_major_formatter(formatter)

    fig = plt.gcf()
    plt.show()
    
    # Make an output filename out of the filename and the variable itself.
    base = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.O3.png".format(os.path.basename(__file__))
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'TES-Aura_L2-O3-Nadir_r0000011015_F05_07.he5'

    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
