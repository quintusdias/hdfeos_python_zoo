"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LAADS NPP VIIRS
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python NPP_D16BRDF3_L3D_A2012241_h20v03_C1_03001_2012258151353.py

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
    # files, otherwise look in the current directory.
    FILE_NAME = 'NPP_D16BRDF3_L3D.A2012241.h20v03.C1_03001.2012258151353.hdf'
    LAT_FILE_NAME = 'lat_NPP_D16BRDF3_L3D.A2012241.h20v03.C1_03001.2012258151353.output'
    LON_FILE_NAME = 'lon_NPP_D16BRDF3_L3D.A2012241.h20v03.C1_03001.2012258151353.output'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)
        LAT_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                     LAT_FILE_NAME)
        LON_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                     LON_FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Albedo_BSA_Band1'

    if USE_GDAL:

        import gdal

        # Read dataset.
        GRID_NAME = 'NPP_Grid_BRDF'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)

        # Read parameters for constructing the grid.
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)

        # Construct the grid.
        x = np.linspace(x0, x0 + xinc*nx, nx)
        y = np.linspace(y0, y0 + yinc*ny, ny)
        xv, yv = np.meshgrid(x, y)

        # In basemap, the sinusoidal projection is global, so we won't use it.
        # Instead we'll convert the grid back to lat/lons.
        sinu = pyproj.Proj("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
        wgs84 = pyproj.Proj("+init=EPSG:4326")
        lon, lat = pyproj.transform(sinu, wgs84, xv, yv)

        # Read fill value, valid range, scale factor, add_offset attributes.
        meta = gdset.GetMetadata()

        # Apply the scale factor, valid range, fill value because GDAL does not
        # do this.  Also, GDAL reads the attributes as character values, so we
        # have to properly convert them.
        _FillValue = float(meta['_FillValue'])
        valid_range = [float(x) for x in meta['valid_range'].split(', ')]
        scale_factor = float(meta['scale_factor'])
        add_offset = float(meta['add_offset'])
        units = meta['units']
        long_name = meta['long_name']

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.double)

        # Read geolocation dataset from HDF-EOS2 dumper output.
        lat = np.genfromtxt(LAT_FILE_NAME, delimiter=',', usecols=[0])
        lat = lat.reshape(data.shape)

        lon = np.genfromtxt(LON_FILE_NAME, delimiter=',', usecols=[0])
        lon = lon.reshape(data.shape)

        # Read attributes
        attrs = data2D.attributes(full=1)
        long_name = attrs["long_name"][0]
        valid_range = attrs["valid_range"][0]
        add_offset = attrs["add_offset"][0]
        _FillValue = attrs["_FillValue"][0]
        scale_factor = attrs["scale_factor"][0]
        units = attrs["units"][0]

    invalid = np.logical_or(data < valid_range[0], data > valid_range[1])
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = data * scale_factor + add_offset
    data = np.ma.masked_array(data, np.isnan(data))

    m = Basemap(projection='cyl', resolution='i',
                lon_0=-10,
                llcrnrlat=45, urcrnrlat=65,
                llcrnrlon=25, urcrnrlon=65)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(45, 61, 5), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(25, 56, 10), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data, latlon=True)

    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    plt.show(block=False)
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
