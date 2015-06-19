"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an LP DAAC MCD43A3
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MCD43A3_A2013305_h12v11_005_2013322102420.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

USE_GDAL = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'MCD43A3.A2013305.h12v11.005.2013322102420.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Albedo_BSA_Band1'

    if USE_GDAL:

        import gdal

        GRID_NAME = 'MOD_Grid_BRDF'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)

        # Construct the grid.
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
        x = np.linspace(x0, x0 + xinc*nx, nx)
        y = np.linspace(y0, y0 + yinc*ny, ny)
        xv, yv = np.meshgrid(x, y)

        # In basemap, the sinusoidal projection is global, so we won't use it.
        # Instead we'll convert the grid back to lat/lons.
        sinu = pyproj.Proj("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
        wgs84 = pyproj.Proj("+init=EPSG:4326")
        lon, lat = pyproj.transform(sinu, wgs84, xv, yv)

        # Read the attributes.
        meta = gdset.GetMetadata()
        long_name = meta['long_name']
        units = meta['units']
        _FillValue = np.float(meta['_FillValue'])
        scale_factor = np.float(meta['scale_factor'])
        add_offset = np.float(meta['add_offset'])
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')]

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.double)

        # Read geolocation dataset from HDF-EOS2 dumper output.
        GEO_FILE_NAME = 'lat_MCD43A3.A2013305.h12v11.005.2013322102420.output'
        GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                     GEO_FILE_NAME)
        lat = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        lat = lat.reshape(data.shape)

        GEO_FILE_NAME = 'lon_MCD43A3.A2013305.h12v11.005.2013322102420.output'
        GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                     GEO_FILE_NAME)
        lon = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        lon = lon.reshape(data.shape)

        # Read attributes.
        attrs = data2D.attributes(full=1)
        long_name = attrs["long_name"][0]
        valid_range = attrs["valid_range"][0]
        add_offset = attrs["add_offset"][0]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["scale_factor"][0]
        units = attrs["units"][0]

    invalid = data == _FillValue
    invalid = np.logical_or(invalid, data < valid_range[0])
    invalid = np.logical_or(invalid, data > valid_range[1])
    data[invalid] = np.nan
    data = scale_factor * (data - add_offset)
    data = np.ma.masked_array(data, np.isnan(data))

    m = Basemap(projection='cyl', resolution='l',
                lon_0=-10,
                llcrnrlat=-32.5, urcrnrlat=-17.5,
                llcrnrlon=-72.5, urcrnrlon=-52.5)
    m.drawcoastlines(linewidth=1.0)
    m.drawparallels(np.arange(-30, -10, 5), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-70, -50, 5), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data)
    # Subset if you want to speed up processing.
    # m.pcolormesh(lon[::2], lat[::2], data[::2])

    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}\n'.format(basename, long_name), fontsize=11)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}\n'.format(basename, long_name), fontsize=11)

    fig = plt.gcf()
    # plt.show()

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
