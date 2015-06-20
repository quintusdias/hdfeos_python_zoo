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

    python AMSR_E_L3_SI_12km_SH_36H_DAY.py

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
    FILE_NAME = 'AMSR_E_L3_SeaIce12km_V11_20050118.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'SI_12km_SH_36H_DAY'

    if USE_GDAL:

        import gdal

        GRID_NAME = 'SpPolarGrid12km'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)

        # Read projection parameters from global attribute.
        meta = gdset.GetMetadata()
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)
        data = hdf.select(DATAFIELD_NAME)[:].astype(np.float64)

        # There are multiple grids in this file.
        # Thus, simple regular expression search for
        # UpperLeft/LowerRightPoint from StructMetadata.0 won't work.
        #
        # Use HDFView and look for the following parameters:
        #
        # GROUP=GRID_2
        #  GridName="SpPolarGrid06km"
        #    XDim=1264
        #    YDim=1328
        #    UpperLeftPointMtrs=(-3950000.000000,4350000.000000)
        #    LowerRightMtrs=(3950000.000000,-3950000.000000)
        ny, nx = data.shape
        x1 = 3950000
        x0 = -3950000
        y0 = 4350000
        y1 = -3950000
        xinc = (x1 - x0) / nx
        yinc = (y1 - y0) / ny

    # Apply the attributes information.
    # Ref:  http://nsidc.org/data/docs/daac/ae_si12_12km_seaice/data.html
    data[data == 0] = np.nan
    data *= 0.1
    data = np.ma.masked_array(data, np.isnan(data))

    x = np.linspace(x0, x0 + xinc*nx, nx)
    y = np.linspace(y0, y0 + yinc*ny, ny)
    xv, yv = np.meshgrid(x, y)
    args = ["+proj=stere",
            "+lat_0=-90",
            "+lon_0=0",
            "+lat_ts=-70",
            "+k=1",
            "+es=0.006693883",
            "+a=6378273",
            "+x_0=0",
            "+y_0=0",
            "+ellps=WGS84",
            "+datum=WGS84"]
    pstereo = pyproj.Proj(' '.join(args))
    wgs84 = pyproj.Proj("+init=EPSG:4326")
    lon, lat = pyproj.transform(pstereo, wgs84, xv, yv)

    units = 'K'
    long_name = DATAFIELD_NAME

    m = Basemap(projection='spstere', resolution='l', boundinglat=-45, lon_0=0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-80, 0, 20), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 30), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data, latlon=True)

    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    plt.show(block=False)
    pngfile = "{0}.1.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
