"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LAADS
MYD swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MYD02HKM_A2010031_0035_005_2010031183706_EV_500_RefSB_lvl0.py

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
    FILE_NAME = 'MYD02HKM.A2010031.0035.005.2010031183706.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'EV_500_RefSB'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)

        # The following method doesn't rely on HDF-EOS2 dumper outputs and
        # plots at the half resolution form the original.
        #
        # Just read the first level, and subset the data to match the lat/lon
        # resolution.
        data = nc.variables[DATAFIELD_NAME][0, ::2, ::2].astype(np.float64)

        # Retrieve geo-location datasets.
        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]

        data = nc.variables[DATAFIELD_NAME][0, :, :].astype(np.float64)

        # Read geolocation dataset from HDF-EOS2 dumper output.
        # GEO_FILE_NAME = 'lat_MYD02HKM.A2010031.0035.005.2010031183706.output'
        # GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
        #                              GEO_FILE_NAME)
        # latitude = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        # latitude = latitude.reshape(data.shape)

        # GEO_FILE_NAME = 'lon_MYD02HKM.A2010031.0035.005.2010031183706.output'
        # GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
        #                              GEO_FILE_NAME)
        # longitude = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        # longitude = longitude.reshape(data.shape)

        units = nc.variables[DATAFIELD_NAME].reflectance_units
        long_name = nc.variables[DATAFIELD_NAME].long_name

        # The scale and offset attributes do not have standard names in this
        # case, so we have to apply the scaling equation ourselves.
        # Fill value is already applied, though.
        scale_factor = nc.variables[DATAFIELD_NAME].reflectance_scales[0]
        add_offset = nc.variables[DATAFIELD_NAME].reflectance_offsets[1]
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

        # The following method doesn't rely on HDF-EOS2 dumper outputs and
        # plots at the half resolution form the original.
        #
        # Just read the first level, and subset the data to match the lat/lon
        # resolution.
        data = data3D[0, ::2, ::2].astype(np.double)
        lat = hdf.select('Latitude')
        latitude = lat[:]
        lon = hdf.select('Longitude')
        longitude = lon[:]

        # data = data3D[0, :, :].astype(np.double)

        # Read geolocation dataset from HDF-EOS2 dumper output.
        # GEO_FILE_NAME = 'lat_MYD02HKM.A2010031.0035.005.2010031183706.output'
        # GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
        #                              GEO_FILE_NAME)
        # latitude = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        # latitude = latitude.reshape(data.shape)

        # GEO_FILE_NAME = 'lon_MYD02HKM.A2010031.0035.005.2010031183706.output'
        # GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
        #                              GEO_FILE_NAME)
        # longitude = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        # longitude = longitude.reshape(data.shape)

        # Retrieve attributes.
        attrs = data3D.attributes(full=1)
        long_name = attrs["long_name"][0]
        add_offset = attrs["reflectance_offsets"][0][1]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["reflectance_scales"][0][1]
        valid_min = attrs["valid_range"][0][0]
        valid_max = attrs["valid_range"][0][1]
        units = attrs["reflectance_units"][0]

        # Retrieve dimension name.
        dim = data3D.dim(0)
        dimname = dim.info()[0]

    invalid = np.logical_or(data > valid_max, data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = scale_factor * (data - add_offset)
    data = np.ma.masked_array(data, np.isnan(data))

    # Subsample data. Otherwise, you will notice that some regions are not
    # plotted on Mac OS X. I suspect that it is a Python's map plot
    # or system's limilt.
    data = data[::2, ::2]
    latitude = latitude[::2, ::2]
    longitude = longitude[::2, ::2]

    # Use a hemispherical projection for the southern hemisphere since the
    # swath is over Antarctica.
    m = Basemap(projection='splaea', resolution='h',
                boundinglat=-65, lon_0=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, -50, 10), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    title = '{0}\n{1}\nat {2}=0'
    title = title.format(basename,
                         'Reflectance derived from ' + long_name, dimname)
    plt.title(title, fontsize=10)
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
