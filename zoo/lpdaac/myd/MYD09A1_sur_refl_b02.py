"""
This example code illustrates how to access and visualize an LP_DAAC MYD
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
"""

import os
import re

import gdal
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

USE_NETCDF = True

def run(FILE_NAME):
    
    DATAFIELD_NAME = 'sur_refl_b02'

    if USE_NETCDF:

        from netCDF4 import Dataset

        # The scaling equation isn't what netcdf4 expects, so turn it off.
        nc = Dataset(FILE_NAME)
        ncvar = nc.variables[DATAFIELD_NAME]
        ncvar.set_auto_maskandscale(False)
        data = ncvar[:].astype(np.float64)

        # Get any needed attributes.
        scale = ncvar.scale_factor
        offset = ncvar.add_offset
        fillvalue = ncvar._FillValue
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
    
    else:
        # Gdal
        import gdal

        GRID_NAME = 'MOD_Grid_500m_Surface_Reflectance'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)
    
        # Get any needed attributes.
        meta = gdset.GetMetadata()
        scale = np.float(meta['scale_factor'])
        offset = np.float(meta['add_offset'])
        fillvalue = np.float(meta['_FillValue'])
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')]
        units = meta['units']
        long_name = meta['long_name']
    
        # Construct the grid.
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
        x = np.linspace(x0, x0 + xinc*nx, nx)
        y = np.linspace(y0, y0 + yinc*ny, ny)
        xv, yv = np.meshgrid(x, y)

        del gdset

    # Apply the attributes to the data.
    invalid = np.logical_or(data < valid_range[0], data > valid_range[1])
    invalid = np.logical_or(invalid, data == fillvalue)
    data[invalid] = np.nan
    data = (data - offset) * scale
    data = np.ma.masked_array(data, np.isnan(data))
    
    # In basemap, the sinusoidal projection is global, so we won't use it.
    # Instead we'll convert the grid back to lat/lons.
    sinu = pyproj.Proj("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
    wgs84 = pyproj.Proj("+init=EPSG:4326") 
    lon, lat= pyproj.transform(sinu, wgs84, xv, yv)

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=7.5, urcrnrlat=22.5,
                llcrnrlon=-162.5, urcrnrlon = -137.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(5, 25, 5), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-170, -130, 10), labels=[0, 0, 0, 1])
    m.pcolormesh(lon[::2], lat[::2], data[::2], latlon=True)
    m.colorbar()
    title = "{0}".format(long_name.replace('_', ' '))
    plt.title(title)

    fig = plt.gcf()
    plt.show()

    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MYD09A1.A2007273.h03v07.005.2007285103507.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
