"""
This example code illustrates how to access and visualize an LP DAAC MCD
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MCD43C1_Black_Sky_Albedo_Num_Albedo_Bands1.py

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
    GRID_NAME = 'MOD_CMG_BRDF_0.05Deg'
    DATAFIELD_NAME = 'Black_Sky_Albedo'
    
    gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                     GRID_NAME,
                                                     DATAFIELD_NAME)

    # Just retrieve the 2nd band.
    gdset = gdal.Open(gname)
    band = gdset.GetRasterBand(2)
    data = band.ReadAsArray().astype(np.float64)

    # Apply the attributes.
    meta = gdset.GetMetadata()
    fillvalue = np.float(meta['_FillValue'])
    scale = np.float(meta['scale_factor'])
    offset = np.float(meta['add_offset'])
    units = meta['units']

    data[data == fillvalue] = np.nan
    data = scale * (data - offset)

    data = np.ma.masked_array(data, np.isnan(data))

    # Normally we would use the following code to reconstruct the grid, but
    # the grid metadata is incorrect in this case, specifically the upper left
    # and lower right coordinates of the grid.  We'll construct the grid
    # manually, taking into account the fact that we're going to subset the
    # data by a factor of 10 (the grid size is 3600 x 7200).
    x = np.linspace(-180, 180, 720)
    y = np.linspace(90, -90, 360)
    lon, lat = np.meshgrid(x, y)

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 90, 45), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data[::10,::10])
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
    hdffile = 'MCD43C1.A2006353.004.2007012185705.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
