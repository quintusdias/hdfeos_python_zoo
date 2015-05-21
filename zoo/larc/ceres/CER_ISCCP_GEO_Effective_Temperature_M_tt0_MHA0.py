"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LaRC CERES file in
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python CER_ISCCP_GEO_Effective_Temperature_M_tt0_MHA0.py

The netCDF file must either be in your current working directory
or in a directory specified by the environment variable HDFEOS_ZOO_DIR.
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
    FILE_NAME = 'CER_ISCCP-D2like-GEO_Composite_Beta1_023031.200510.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Effective Temperature - M'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)

        # Subset the data to match the size of the swath geolocation fields.
        # Turn off autoscaling, we'll handle that ourselves due to presence of
        # a valid range.
        var = nc.variables[DATAFIELD_NAME]
        data = var[0, 0, :, :].astype(np.float64)

        # Read the geolocation.
        longitude = nc.variables['Longitude - MH'][0, :, :].astype(np.float64)
        colatitude = nc.variables['Colatitude - MH'][0, :, :].astype(np.float64)

        # Read attributes.
        fillvalue = var._FillValue
        units = var.units

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data4D = hdf.select(DATAFIELD_NAME)
        data = data4D[0, 0, :, :].astype(np.float64)

        # Read geolocation dataset.
        lon3D = hdf.select('Longitude - MH')
        longitude = lon3D[0, :, :].astype(np.float64)
        lat3D = hdf.select('Colatitude - MH')
        colatitude = lat3D[0, :, :].astype(np.float64)

        # Read attributes.
        attrs = data4D.attributes(full=1)
        fillvalue = attrs["_FillValue"][0]
        units = attrs["units"][0]

    # Apply the attributes.
    data[data == fillvalue] = np.nan
    datam = np.ma.masked_array(data, mask=np.isnan(data))

    latitude = 90 - colatitude

    # The data is global, so render in a global projection.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 90, 45))
    m.drawmeridians(np.arange(-180, 180, 45),
                    labels=[True, False, False, True])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    title = '{0}\n{1}'.format(basename,
                              'Monthly Mean Effective Temperature of Cumulus')
    plt.title(title)
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
