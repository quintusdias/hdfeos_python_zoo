"""
This example code illustrates how to access and visualize a LaRC MOPITT grid
HDF-EOS2 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOP03_CO_Profiles_Day_horizontal_lvl111.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):
    
    nc = Dataset(FILE_NAME)

    DATAFIELD_NAME = 'CO Profiles Day'
    data = nc.variables[DATAFIELD_NAME][111, :, :].astype(np.float64)
    lat = nc.variables['Latitude'][:]
    lon = nc.variables['Longitude'][:]
    pres = nc.variables['Pressure Grid'][:]

    # Replace the fill value with NaN
    data[data == -9999] = np.nan
    data = np.ma.masked_array(data, np.isnan(data))
    
    # Contour the data on a grid of longitude vs. pressure
    longitude, pressure = np.meshgrid(lon, pres)
    plt.contourf(longitude, pressure, data.T)
    plt.title("{0} (PPVB)".format(DATAFIELD_NAME))
    plt.colorbar()

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME.replace(' ', '_'))
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOP03-20000303-L3V1.0.1.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)

