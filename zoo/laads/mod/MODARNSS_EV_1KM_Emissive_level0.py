"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LAADS MODIS swath
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MODARNSS_EV_1KM_Emissive_level0.py

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
    FILE_NAME = 'MODARNSS.Abracos_Hill.A2000080.1515.005.2007164153544.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'EV_1KM_Emissive'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        var = nc.variables[DATAFIELD_NAME]

        # Have to be very careful of the scaling equation here.
        # We'll turn autoscaling off in order to correctly scale the data.
        var.set_auto_maskandscale(False)
        data = var[0, :, :].astype(np.double)

        # Retrieve the geolocation data.
        longitude = nc.variables['Longitude'][:]
        latitude = nc.variables['Latitude'][:]

        # Retrieve attributes.
        scale_factor = var.radiance_scales[0]
        add_offset = var.radiance_offsets[0]
        _FillValue = var._FillValue
        valid_min = var.valid_range[0]
        valid_max = var.valid_range[1]
        long_name = var.long_name
        units = var.radiance_units

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data3D = hdf.select(DATAFIELD_NAME)
        data = data3D[0, :, :].astype(np.double)

        # Read geolocation dataset.
        lat = hdf.select('Latitude')
        latitude = lat[:]
        lon = hdf.select('Longitude')
        longitude = lon[:]

        # Retrieve attributes.
        attrs = data3D.attributes(full=1)
        long_name = attrs["long_name"][0]
        add_offset = attrs["radiance_offsets"][0][0]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["radiance_scales"][0][0]
        valid_min = attrs["valid_range"][0][0]
        valid_max = attrs["valid_range"][0][1]
        units = attrs["radiance_units"][0]

    invalid = np.logical_or(data > valid_max,
                            data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = (data - add_offset) * scale_factor
    data = np.ma.masked_array(data, np.isnan(data))

    # Render the plot in a cylindrical projection.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-12, urcrnrlat=-9,
                llcrnrlon=-64, urcrnrlon=-61)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-12., -8., 1.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-64, -60., 1), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}\n'.format(basename,
                                  'Radiance derived from ' + long_name),
              fontsize=11)
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
