"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LaRC MOPITT grid
HDF-EOS2 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOP03_CO_Profiles_Day_lvl1.py

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
    FILE_NAME = 'MOP03-20000303-L3V1.0.1.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'CO Profiles Day'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)

        data = nc.variables[DATAFIELD_NAME][:, :, 1].astype(np.float64)
        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]
        pressure = nc.variables['Pressure Grid'][:]

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        data = hdf.select(DATAFIELD_NAME)[:, :, 1].astype(np.float64)
        latitude = hdf.select('Latitude')[:]
        longitude = hdf.select('Longitude')[:]
        pressure = hdf.select('Pressure Grid')[:]

    # Replace the fill value with NaN
    data[data == -9999] = np.nan
    data = np.ma.masked_array(data, np.isnan(data))

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 45), labels=[True, False, False, True])
    m.drawmeridians(np.arange(-180, 180, 45),
                    labels=[True, False, False, True])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    cb = m.colorbar()
    cb.set_label('ppbv')

    basename = os.path.basename(FILE_NAME)
    title = '{0}\n{1} at Pressure={2} hPa'.format(basename,
                                                  DATAFIELD_NAME,
                                                  pressure[1])
    plt.title(title)

    fig = plt.gcf()
    plt.show(block=False)
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
