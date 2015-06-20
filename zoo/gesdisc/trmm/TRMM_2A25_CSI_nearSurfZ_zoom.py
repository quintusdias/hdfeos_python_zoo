"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC TRMM 2A25
 version 6 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python TRMM_2A25_CSI_nearSurfZ_zoom.py

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
    FILE_NAME = '2A25_CSI.990804.9692.KORA.6.HDF'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'nearSurfZ'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        data = nc.variables[DATAFIELD_NAME][:].astype(np.float64)

        # Retrieve the geolocation data.
        latitude = nc.variables['geolocation'][:, :, 0]
        longitude = nc.variables['geolocation'][:, :, 1]

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)
        ds = hdf.select(DATAFIELD_NAME)
        data = ds[:].astype(np.double)

        # Retrieve the geolocation data.
        geo = hdf.select('geolocation')
        latitude = geo[:, :, 0]
        longitude = geo[:, :, 1]

    # There's no fill value set, but 0.0 is considered the fill value.
    data[data == 0.0] = np.nan
    datam = np.ma.masked_array(data, np.isnan(data))

    # Draw an equidistant cylindrical projection using the high resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='h',
                llcrnrlat=30, urcrnrlat=36,
                llcrnrlon=123, urcrnrlon=135)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(30, 37), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(123, 135, 2), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label('Units: none')

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, DATAFIELD_NAME))
    fig = plt.gcf()
    plt.show(block=False)

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
