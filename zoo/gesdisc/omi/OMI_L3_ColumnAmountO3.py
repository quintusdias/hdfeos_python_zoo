"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC OMI L3 file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python OMI_L3_ColumnAmountO3.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

FILE_NAME = 'OMI-Aura_L3-OMTO3e_2005m1214_v002-2006m0929t143855.he5'

# Can do this using either netCDF4 or h5py.
USE_NETCDF4 = True


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'OMI-Aura_L3-OMTO3e_2005m1214_v002-2006m0929t143855.he5'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    if USE_NETCDF4:

        from netCDF4 import Dataset

        DATAFIELD_NAME = 'ColumnAmountO3'
        nc = Dataset(FILE_NAME)
        grp = nc.groups['HDFEOS'].groups['GRIDS']
        grp = grp.groups['OMI Column Amount O3'].groups['Data Fields']
        var = grp.variables[DATAFIELD_NAME]
        data = var[:]

        # Get attributes needed for the plot.
        title = var.Title
        units = var.Units
        nc.close()

    else:

        import h5py

        DATAFIELD_NAME = '/HDFEOS/GRIDS/OMI Column Amount O3/Data Fields/ColumnAmountO3'
        with h5py.File(FILE_NAME, mode='r') as f:

            dset = f[DATAFIELD_NAME]
            data = dset[:]

            # Have to manually create a masked array due to the fill value.
            # No need to scale the data, as the scale factor and add offset are
            # 1.0 and 0.0 respectively.
            data[data == dset.fillvalue] = np.nan
            data = np.ma.masked_where(np.isnan(data), data)

            # Get attributes needed for the plot.
            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            title = dset.attrs['Title'].decode()
            units = dset.attrs['Units'].decode()

    # There is no geolocation data, so construct it ourselves.
    longitude = np.arange(0., 1440.0) * 0.25 - 180 + 0.125
    latitude = np.arange(0., 720.0) * 0.25 - 90 + 0.125

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
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, title))

    fig = plt.gcf()
    # plt.show()

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
