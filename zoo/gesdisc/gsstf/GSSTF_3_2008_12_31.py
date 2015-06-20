"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC MEaSUREs
GSSTF HDF-EOS5 grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python GSSTF_3_2008_12_31.py

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
    FILE_NAME = 'GSSTF.3.2008.12.31.he5'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        grp = nc.groups['HDFEOS'].groups['GRIDS'].groups['SET1']
        data_var = grp.groups['Data Fields'].variables['E']
        data = data_var[:]

        data_longname = data_var.long_name
        data_units = data_var.units

    else:

        import h5py

        with h5py.File(FILE_NAME, mode='r') as f:

            dset_var = f['/HDFEOS/GRIDS/SET1/Data Fields/E']
            data = dset_var[:]

            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            data_units = dset_var.attrs['units'].decode()
            data_longname = dset_var.attrs['long_name'].decode()
            fv = dset_var.attrs['_FillValue'][0]

            # We have to apply the fill value ourselves.
            data[data == fv] = np.nan
            data = np.ma.masked_array(data, np.isnan(data))

    # The projection is GEO, so we can construct the lat/lon arrays ourselves.
    scaleX = 360.0 / data.shape[1]
    scaleY = 180.0 / data.shape[0]
    longitude = np.arange(data.shape[1]) * scaleX - 180 + scaleX/2
    latitude = np.arange(data.shape[0]) * scaleY - 90 + scaleY/2

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(data_units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, data_longname))

    fig = plt.gcf()
    plt.show(block=False)

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
