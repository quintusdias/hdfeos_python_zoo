"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an NSIDC AMSR grid
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AMSR_E_L3_DL_A_TB36_5H_Res_1.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_GDAL = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'AMSR_E_L3_DailyLand_V06_20050118.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'A_TB36.5H (Res 1)'

    if USE_GDAL:

        import gdal
        import mpl_toolkits.basemap.pyproj as pyproj

        GRID_NAME = 'Ascending_Land_Grid'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)

        meta = gdset.GetMetadata()
        _FillValue = float(meta['_FillValue'])

        # Construct the grid.
        # Reproject out of the global GCTP CEA into lat/lon.
        # Ref:  http://nsidc.org/data/atlas/epsg_3410.html
        meta = gdset.GetMetadata()
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
        x = np.linspace(x0, x0 + xinc*nx, nx)
        y = np.linspace(y0, y0 + yinc*ny, ny)
        xv, yv = np.meshgrid(x, y)
        args = ["+proj=cea",
                "+lat_0=0",
                "+lon_0=0",
                "+lat_ts=30",
                "+a=6371228",
                "+units=m"]
        pstereo = pyproj.Proj(' '.join(args))
        wgs84 = pyproj.Proj("+init=EPSG:4326")
        lon, lat = pyproj.transform(pstereo, wgs84, xv, yv)
        del gdset

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.float64)

        # Read geolocation dataset from HDF-EOS2 dumper output.
        base = 'AMSR_E_L3_DailyLand_V06_20050118_Ascending_Land_Grid.output'
        GEO_FILE_NAME = 'lat_' + base
        try:
            GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                         GEO_FILE_NAME)
        except KeyError:
            pass

        lat = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        GEO_FILE_NAME = 'lon_' + base
        try:
            GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                         GEO_FILE_NAME)
        except KeyError:
            pass
        lon = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])

        _FillValue = data2D.attributes(full=1)["_FillValue"][0]

    # Apply the attributes information.
    # Ref:  http://nsidc.org/data/docs/daac/ae_land3_l3_soil_moisture/data.html
    data[data == _FillValue] = np.nan
    data *= 0.1
    data = np.ma.masked_array(data, np.isnan(data))
    long_name = DATAFIELD_NAME
    units = 'Kelvin'

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, llcrnrlon=-180, urcrnrlat=90, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
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
