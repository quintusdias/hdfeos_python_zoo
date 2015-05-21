"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LAADS MYD swath
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MYDARNSS_EV_1KM_Emissive_lvl9.py

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
    # files, otherwise look in the current directory.
    FILE_NAME = 'MYDARNSS.Barrow.A2002184.2200.005.2007051063709.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'EV_1KM_Emissive'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)

        # Read 3D dataset and subset it.
        data = nc.variables[DATAFIELD_NAME][9, :, :].astype(np.float64)

        # Read geo-location dataset.
        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]

        # Read attributes.
        units = nc.variables[DATAFIELD_NAME].radiance_units
        long_name = nc.variables[DATAFIELD_NAME].long_name

        # The scale and offset attributes do not have standard names in this
        # case, so we have to apply the scaling equation ourselves.
        scale_factor = nc.variables[DATAFIELD_NAME].radiance_scales[9]
        add_offset = nc.variables[DATAFIELD_NAME].radiance_offsets[9]
        valid_range = nc.variables[DATAFIELD_NAME].valid_range
        _FillValue = nc.variables[DATAFIELD_NAME]._FillValue
        valid_min = valid_range[0]
        valid_max = valid_range[1]

        # Retrieve dimension name.
        dimname = nc.variables[DATAFIELD_NAME].dimensions[0]

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data3D = hdf.select(DATAFIELD_NAME)
        data = data3D[9, :, :].astype(np.double)

        # Read geolocation dataset.
        lat = hdf.select('Latitude')
        latitude = lat[:]
        lon = hdf.select('Longitude')
        longitude = lon[:]

        # Retrieve attributes.
        attrs = data3D.attributes(full=1)
        long_name = attrs["long_name"][0]
        add_offset = attrs["radiance_offsets"][0][9]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["radiance_scales"][0][9]
        valid_min = attrs["valid_range"][0][0]
        valid_max = attrs["valid_range"][0][1]
        units = attrs["radiance_units"][0]

        # Retrieve dimension name.
        dim = data3D.dim(0)
        dimname = dim.info()[0]

    invalid = np.logical_or(data > valid_max,
                            data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = scale_factor * (data - add_offset)
    data = np.ma.masked_array(data, np.isnan(data))

    m = Basemap(projection='laea', resolution='i',
                lat_ts=71.25, lat_0=71.25, lon_0=-156.5,
                width=100000, height=100000)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(70, 72.1, 0.5), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-158, -154.9, 0.5), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    title = '{0}\n{1}\nat {2}=9'
    title = title.format(basename,
                         'Radiance derived from ' + long_name, dimname)
    plt.title(title, fontsize=11)

    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
