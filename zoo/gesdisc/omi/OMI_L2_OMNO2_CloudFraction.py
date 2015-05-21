"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC OMI L2 file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python OMI_L2_OMNO2_CloudFraction.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = True


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'OMI-Aura_L2-OMNO2_2008m0720t2016-o21357_v003-2008m0721t101450.he5'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'CloudFraction'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        grp = nc.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2']
        var = grp.groups['Data Fields'].variables[DATAFIELD_NAME]

        # netCDF4 doesn't quite handle the scaling correctly in this case since
        # the attributes are "ScaleFactor" and "AddOffset" instead of
        # "scale_factor" and "add_offset".  We'll do the scaling and conversion
        # to a masked array ourselves.
        var.set_auto_maskandscale(False)

        data = var[:].astype(np.float64)

        # Retrieve any attributes that may be needed later.
        scale = var.ScaleFactor
        offset = var.Offset
        title = var.Title
        missing_value = var.MissingValue
        fill_value = var._FillValue
        units = var.Units

        # Retrieve the geolocation data.
        latitude = grp.groups['Geolocation Fields'].variables['Latitude'][:]
        longitude = grp.groups['Geolocation Fields'].variables['Longitude'][:]

    else:

        import h5py

        path = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/'
        DATAFIELD_NAME = path + 'CloudFraction'

        with h5py.File(FILE_NAME, mode='r') as f:

            dset = f[DATAFIELD_NAME]
            data = dset[:].astype(np.float64)

            # Retrieve any attributes that may be needed later.
            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            scale = f[DATAFIELD_NAME].attrs['ScaleFactor']
            offset = f[DATAFIELD_NAME].attrs['Offset']
            missing_value = f[DATAFIELD_NAME].attrs['MissingValue']
            fill_value = f[DATAFIELD_NAME].attrs['_FillValue']
            title = f[DATAFIELD_NAME].attrs['Title'].decode()
            units = f[DATAFIELD_NAME].attrs['Units'].decode()

            # Retrieve the geolocation data.
            path = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/'
            latitude = f[path + 'Latitude'][:]
            longitude = f[path + 'Longitude'][:]

    data[data == missing_value] = np.nan
    data[data == fill_value] = np.nan
    data = scale * (data - offset)
    datam = np.ma.masked_where(np.isnan(data), data)

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, title), fontsize=12)

    fig = plt.gcf()
    # plt.show()

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
