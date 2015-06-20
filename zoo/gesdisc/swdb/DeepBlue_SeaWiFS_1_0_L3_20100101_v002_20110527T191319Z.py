"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC MEaSUREs
SeaWiFS L3 grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python DeepBlue_SeaWiFS_1_0_L3_20100101_v002_20110527T191319Z.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

# Can do this using either netCDF4 or h5py.
USE_NETCDF4 = True


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'DeepBlue-SeaWiFS-1.0_L3_20100101_v002-20110527T191319Z.h5'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'aerosol_optical_thickness_550_ocean'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        data = nc.variables[DATAFIELD_NAME][:]
        latitude = nc.variables['latitude'][:]
        longitude = nc.variables['longitude'][:]

        # Get attributes needed for the plot.
        long_name = nc.variables[DATAFIELD_NAME].long_name

    else:

        import h5py

        with h5py.File(FILE_NAME, mode='r') as f:

            data = f[DATAFIELD_NAME][:]
            latitude = f['latitude'][:]
            longitude = f['longitude'][:]

            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            long_name = f[DATAFIELD_NAME].attrs['long_name'].decode()

    data[data == -999] = np.nan
    data = np.ma.masked_array(data, np.isnan(data))

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(-180, 180., 45.))
    m.pcolormesh(longitude, latitude, data, vmin=0, vmax=1, latlon=True)
    m.colorbar()

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, long_name))

    fig = plt.gcf()
    plt.show(block=False)

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
