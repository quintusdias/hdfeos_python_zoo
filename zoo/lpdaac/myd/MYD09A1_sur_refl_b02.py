"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an LP_DAAC MYD09A1
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MYD09A1_sur_refl_b02.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

In order for the netCDF code path to work, the netcdf library must be compiled
with HDF4 support.  Please see the README for details.
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

USE_NETCDF4 = False
USE_GDAL = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'MYD09A1.A2007273.h03v07.005.2007285103507.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'sur_refl_b02'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        # The scaling equation isn't what netcdf4 expects, so turn it off.
        nc = Dataset(FILE_NAME)
        ncvar = nc.variables[DATAFIELD_NAME]
        ncvar.set_auto_maskandscale(False)
        data = ncvar[:].astype(np.float64)

        # Get any needed attributes.
        scale_factor = ncvar.scale_factor
        add_offset = ncvar.add_offset
        _FillValue = ncvar._FillValue
        valid_range = ncvar.valid_range
        units = ncvar.units
        long_name = ncvar.long_name

        # Construct the grid.  The needed information is in a global attribute
        # called 'StructMetadata.0'.  Use regular expressions to tease out the
        # extents of the grid.
        gridmeta = getattr(nc, 'StructMetadata.0')
        ul_regex = re.compile(r'''UpperLeftPointMtrs=\(
                                  (?P<upper_left_x>[+-]?\d+\.\d+)
                                  ,
                                  (?P<upper_left_y>[+-]?\d+\.\d+)
                                  \)''', re.VERBOSE)
        match = ul_regex.search(gridmeta)
        x0 = np.float(match.group('upper_left_x'))
        y0 = np.float(match.group('upper_left_y'))

        lr_regex = re.compile(r'''LowerRightMtrs=\(
                                  (?P<lower_right_x>[+-]?\d+\.\d+)
                                  ,
                                  (?P<lower_right_y>[+-]?\d+\.\d+)
                                  \)''', re.VERBOSE)
        match = lr_regex.search(gridmeta)
        x1 = np.float(match.group('lower_right_x'))
        y1 = np.float(match.group('lower_right_y'))

        nx, ny = data.shape
        x = np.linspace(x0, x1, nx)
        y = np.linspace(y0, y1, ny)
        xv, yv = np.meshgrid(x, y)

        # In basemap, the sinusoidal projection is global, so we won't use it.
        # Instead we'll convert the grid back to lat/lons.
        sinu = pyproj.Proj("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
        wgs84 = pyproj.Proj("+init=EPSG:4326")
        lon, lat = pyproj.transform(sinu, wgs84, xv, yv)

    elif USE_GDAL:

        import gdal

        GRID_NAME = 'MOD_Grid_500m_Surface_Reflectance'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)

        # Get any needed attributes.
        meta = gdset.GetMetadata()
        scale_factor = np.float(meta['scale_factor'])
        add_offset = np.float(meta['add_offset'])
        _FillValue = np.float(meta['_FillValue'])
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')]
        units = meta['units']
        long_name = meta['long_name']

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

        del gdset

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.double)

        # Read geolocation dataset from HDF-EOS2 dumper output.
        GEO_FILE_NAME = 'lat_MYD09A1.A2007273.h03v07.005.2007285103507.output'
        GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                     GEO_FILE_NAME)
        lat = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        lat = lat.reshape(data.shape)

        GEO_FILE_NAME = 'lon_MYD09A1.A2007273.h03v07.005.2007285103507.output'
        GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                     GEO_FILE_NAME)
        lon = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        lon = lon.reshape(data.shape)

        # Read attributes.
        attrs = data2D.attributes(full=1)
        long_name = attrs["long_name"][0]
        valid_range = attrs["valid_range"][0]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["scale_factor"][0]
        add_offset = attrs["add_offset"][0]
        units = attrs["units"][0]

    # Apply the attributes to the data.
    invalid = np.logical_or(data < valid_range[0], data > valid_range[1])
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = (data - add_offset) * scale_factor
    data = np.ma.masked_array(data, np.isnan(data))

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=7.5, urcrnrlat=22.5,
                llcrnrlon=-162.5, urcrnrlon=-137.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(5, 25, 5), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-170, -130, 10), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
