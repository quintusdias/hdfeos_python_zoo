"""
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
    GRID_NAME = 'NPP_Grid_BRDF'
    DATAFIELD_NAME = 'Albedo_BSA_Band1'
    
    gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                     GRID_NAME,
                                                     DATAFIELD_NAME)
    gdset = gdal.Open(gname)

    data = gdset.ReadAsArray().astype(np.float64)
    meta = gdset.GetMetadata()

    # Apply the scale factor, valid range, fill value because GDAL does not
    # do this.  Also, GDAL reads the attributes as character values, so we have
    # to properly convert them.
    fill_value = float(meta['_FillValue'])
    data[data == fill_value] = np.nan
    valid_range = [float(x) for x in meta['valid_range'].split(', ')]
    data[data < valid_range[0]] = np.nan
    data[data > valid_range[1]] = np.nan
    scale_factor = float(meta['scale_factor'])
    data = data * scale_factor

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

    m = Basemap(projection='cyl', resolution='i',
                lon_0=-10,
                llcrnrlat=45, urcrnrlat = 65,
                llcrnrlon=25, urcrnrlon = 65)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(45, 61, 5), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(25, 56, 10), labels=[0, 0, 0, 1])

    m.pcolormesh(lon, lat, data, latlon=True)
    m.colorbar()
    
    fig = plt.gcf()
    
    plt.title(DATAFIELD_NAME.replace('_', ' '))
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)

    del gdset


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'NPP_D16BRDF3_L3D.A2012241.h20v03.C1_03001.2012258151353.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
