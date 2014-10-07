"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an LP_DAAC MEaSURES
WELD CONUS Albers grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python CONUS_annual_2012_h01v06_doy007to356_v1_5.py

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
from netCDF4 import Dataset
import numpy as np

USE_GDAL = False
USE_NETCDF = False

def run(FILE_NAME):
    
    DATAFIELD_NAME = 'NDVI_TOA'
    if USE_GDAL:    

    # Scale down the data by a factor of 5 so that low-memory machines
    # can handle it.
    nc = Dataset(FILE_NAME)
    ncvar = nc.variables[DATAFIELD_NAME]
    ncvar.set_auto_maskandscale(False)
    data = ncvar[::5, ::5].astype(np.float64)

    # Get any needed attributes.
    scale = ncvar.scale_factor
    fillvalue = ncvar._FillValue
    valid_range = ncvar.valid_range
    units = ncvar.units

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
    
    ny, nx = data.shape
    x = np.linspace(x0, x1, nx)
    y = np.linspace(y0, y1, ny)
    xv, yv = np.meshgrid(x, y)
    

        # Gdal
        import gdal

        GRID_NAME = 'WELD_GRID'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)

        # Scale down the data by a factor of 5 so that low-memory machines
        # can handle it.
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)[::5, ::5]
    
        # Get any needed attributes.
        meta = gdset.GetMetadata()
        scale = np.float(meta['scale_factor'])
        fillvalue = np.float(meta['_FillValue'])
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')]
        units = meta['units']
    
        # Construct the grid.
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        ny, nx = (gdset.RasterYSize / 5, gdset.RasterXSize / 5)
        x = np.linspace(x0, x0 + xinc*5*nx, nx)
        y = np.linspace(y0, y0 + yinc*5*ny, ny)
        xv, yv = np.meshgrid(x, y)

        del gdset

    else:
        if USE_NETCDF:

            from netCDF4 import Dataset

            # Scale down the data by a factor of 5 so that low-memory machines
            # can handle it.
            nc = Dataset(FILE_NAME)
            ncvar = nc.variables[DATAFIELD_NAME]
            ncvar.set_auto_maskandscale(False)
            data = ncvar[::5, ::5].astype(np.float64)

            # Get any needed attributes.
            scale = ncvar.scale_factor
            fillvalue = ncvar._FillValue
            valid_range = ncvar.valid_range
            units = ncvar.units
            gridmeta = getattr(nc, 'StructMetadata.0')
        else:
            from pyhdf.SD import SD, SDC
            hdf = SD(FILE_NAME, SDC.READ)

            # Read dataset.
            data2D = hdf.select(DATAFIELD_NAME)
            data = data2D[:,:].astype(np.double)

            # Scale down the data by a factor of 6 so that low-memory machines
            # can handle it.
            data = data[::6, ::6]

        
            # Read attributes.
            attrs = data2D.attributes(full=1)
            vra=attrs["valid_range"]
            valid_range = vra[0]
            fva=attrs["_FillValue"]
            fillvalue = fva[0]
            sfa=attrs["scale_factor"]
            scale = sfa[0]        
            ua=attrs["units"]
            units = ua[0]
            fattrs = hdf.attributes(full=1)
            ga = fattrs["StructMetadata.0"]
            gridmeta = ga[0]
        # Construct the grid.  The needed information is in a global attribute
        # called 'StructMetadata.0'.  Use regular expressions to tease out the
        # extents of the grid.  

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
        
        ny, nx = data.shape
        x = np.linspace(x0, x1, nx)
        y = np.linspace(y0, y1, ny)
        xv, yv = np.meshgrid(x, y)

    # Apply the attributes to the data.
    invalid = np.logical_or(data < valid_range[0], data > valid_range[1])
    invalid = np.logical_or(invalid, data == fillvalue)
    data[invalid] = np.nan
    data = data * scale
    data = np.ma.masked_array(data, np.isnan(data))
    
    # Convert the grid back to lat/lon.  The 1st and 2nd standard parallels,
    # the center meridian, and the latitude of projected origin are in the
    # projection parameters contained in the "StructMetadata.0" global
    # attribute.  The following regular expression could have been used to 
    # retrieve them.
    #
    # Ref:  HDF-EOS Library User's Guide for the EOSDIS Evolution and
    #       Development (EED) Contract, Volume 2, Revision 02:  Function
    #       Reference Guide, pages 1-6 through 1-13.
    #
    # aea_regex = re.compile(r'''Projection=GCTP_ALBERS\s+ProjParams=
    #                            \(\d+,\d+
    #                            (?P<stdpr1>[-]?\d+),
    #                            (?P<stdpr2>[-]?\d+),
    #                            (?P<centmer>[-]?\d+),
    #                            (?P<origlat>[-]?\d+)
    #                            ,0{7}\)''', re.VERBOSE)
    #
    aea = pyproj.Proj("+proj=aea +lat_1=29.5 +lat2=45.5 +lon_0=-96 +lat_0=23")
    wgs84 = pyproj.Proj("+init=EPSG:4326") 
    lon, lat= pyproj.transform(aea, wgs84, xv, yv)

    m = Basemap(projection='aea', resolution='i',
                lat_1=29.5, lat_2=45.5, lon_0=-96, lat_0=23,
                llcrnrlat=37.5, urcrnrlat = 42.5,
                llcrnrlon=-127.5, urcrnrlon = -122.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(35, 45, 1), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-130, -120, 1), labels=[0, 0, 0, 1])

    m.pcolormesh(lon, lat, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    long_name = DATAFIELD_NAME
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'CONUS.annual.2012.h01v06.doy007to356.v1.5.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
