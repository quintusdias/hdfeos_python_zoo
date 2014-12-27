"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an NSIDC AMSR_E grid
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AMSR_E_L3_5DaySnow_NH_SWE.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_PYHDFEOS = True
USE_GDAL = False

def run(FILE_NAME):
    
    # Identify the data field.
    DATAFIELD_NAME = 'SWE_NorthernPentad'

    if USE_PYHDFEOS:

        from pyhdfeos import GridFile
        gdf = GridFile(FILE_NAME)
        grid = gdf.grids['Northern Hemisphere']
        data = grid.fields[DATAFIELD_NAME][:].astype(np.float64)
        lat, lon = grid[:]

    elif USE_GDAL:

        import gdal
        import mpl_toolkits.basemap.pyproj as pyproj
        GRID_NAME = 'Northern Hemisphere'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)
        # Construct the grid.
        meta = gdset.GetMetadata()
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
        x = np.linspace(x0, x0 + xinc*nx, nx)
        y = np.linspace(y0, y0 + yinc*ny, ny)
        xv, yv = np.meshgrid(x, y)

        # Reproject the coordinates out of lamaz into lat/lon.
        lamaz = pyproj.Proj("+proj=laea +a=6371228 +lat_0=90 +lon_0=0 +units=m")
        wgs84 = pyproj.Proj("+init=EPSG:4326") 
        lon, lat= pyproj.transform(lamaz, wgs84, xv, yv)
        del gdset

    else:

        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)
        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:,:].astype(np.float64)

	    # Read geolocation dataset from HDF-EOS2 dumper output.
        GEO_FILE_NAME = 'lat_AMSR_E_L3_5DaySnow_V09_20050126.Northern_Hemisphere.output'
        try:
            GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], 
                                         GEO_FILE_NAME)
        except KeyError:
            pass

        lat = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        lat = lat.reshape(data.shape)
        GEO_FILE_NAME = 'lon_AMSR_E_L3_5DaySnow_V09_20050126.Northern_Hemisphere.output'
        try: 
            GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], 
                                         GEO_FILE_NAME)
        except KeyError:
            pass
        lon = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        lon = lon.reshape(data.shape)

    # Filter out invalid range values, multiply by two according to the data
    # spec.
    data[data > 240] = np.nan
    data *= 2
    data = np.ma.masked_array(data, np.isnan(data))
    long_name = 'Northern Hemisphere 5-day Snow Water Equivalent (' + DATAFIELD_NAME + ')'
    units = 'mm'

    # Draw a polar stereographic projection using the low resolution coastline
    # database.
    m = Basemap(projection='npstere', resolution='l',
                boundinglat=25, lon_0 = 0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(0, 91, 20), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180, 30), labels=[0, 0, 0, 1])
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

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'AMSR_E_L3_5DaySnow_V09_20050126.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
