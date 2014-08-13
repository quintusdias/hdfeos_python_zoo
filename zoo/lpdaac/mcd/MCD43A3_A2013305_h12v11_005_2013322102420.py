"""
This example code illustrates how to access and visualize an LP DAAC MODIS
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

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""

import os
import re

import gdal
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

def run(FILE_NAME):
    
    # Identify the data field.
    GRID_NAME = 'MOD_Grid_BRDF'
    DATAFIELD_NAME = 'Albedo_BSA_Band1'
    
    gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                     GRID_NAME,
                                                     DATAFIELD_NAME)
    gdset = gdal.Open(gname)
    data = gdset.ReadAsArray().astype(np.float64)

    # Apply the attributes.
    meta = gdset.GetMetadata()
    units = meta['units']
    fillvalue = np.float(meta['_FillValue'])
    scale = np.float(meta['scale_factor'])
    offset = np.float(meta['add_offset'])
    valid_range = [np.float(x) for x in meta['valid_range'].split(', ')] 

    invalid = data == fillvalue
    invalid = np.logical_or(invalid, data < valid_range[0])
    invalid = np.logical_or(invalid, data > valid_range[1])
    data[invalid] = np.nan
    data = scale * (data - offset)

    data = np.ma.masked_array(data, np.isnan(data))

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
    lon, lat= pyproj.transform(sinu, wgs84, xv, yv)

    m = Basemap(projection='cyl', resolution='l',
                lon_0=-10,
                llcrnrlat=-32.5, urcrnrlat = -17.5,
                llcrnrlon=-72.5, urcrnrlon = -52.5)
    m.drawcoastlines(linewidth=1.0)
    m.drawparallels(np.arange(-30, -10, 5), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-70, -50, 5), labels=[0, 0, 0, 1])
    m.pcolormesh(lon[::2], lat[::2], data[::2])
    m.colorbar()
    title = "{0} ({1})".format(DATAFIELD_NAME.replace('_', ' '), units)
    plt.title(title)

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)

    del gdset


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MCD43A3.A2013305.h12v11.005.2013322102420.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
