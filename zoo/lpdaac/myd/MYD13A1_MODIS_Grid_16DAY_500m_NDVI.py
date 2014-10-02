"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an LP_DAAC MYD13A1
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MYD13A1_MODIS_Grid_16DAY_500m_NDVI.py

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
USE_NETCDF = False

def run(FILE_NAME):
    
    DATAFIELD_NAME = '500m 16 days NDVI'
    if USE_GDAL:
        # Gdal
        import gdal

        GRID_NAME = 'MODIS_Grid_16DAY_500m_VI'
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)

        # Subset the data by a factor of 4 so that low-memory machines can 
        # render it more easily.
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)[::4, ::4]
    
        # Get any needed attributes.
        meta = gdset.GetMetadata()
        scale_factor = np.float(meta['scale_factor'])
        add_offset = np.float(meta['add_offset'])
        _FillValue = np.float(meta['_FillValue'])
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')]
        units = meta['units']
        long_name = meta['long_name']
    
        # Construct the grid, remembering to subset by a factor of 4.
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
        x = np.linspace(x0, x0 + xinc*nx, nx)
        y = np.linspace(y0, y0 + yinc*ny, ny)
        xv, yv = np.meshgrid(x[::4], y[::4])

        del gdset

    else:
        if USE_NETCDF:
            from netCDF4 import Dataset

            # The scaling equation isn't "scale * data + offset", 
            # so turn automatic scaling off.
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
            gridmeta = getattr(nc, 'StructMetadata.0')

        else:
            from pyhdf.SD import SD, SDC
            hdf = SD(FILE_NAME, SDC.READ)

            # Read dataset.
            data2D = hdf.select(DATAFIELD_NAME)
            data = data2D[:,:].astype(np.double)

        
            # Read attributes.
            attrs = data2D.attributes(full=1)
            lna=attrs["long_name"]
            long_name = lna[0]
            vra=attrs["valid_range"]
            valid_range = vra[0]
            aoa=attrs["add_offset"]
            add_offset = aoa[0]
            fva=attrs["_FillValue"]
            _FillValue = fva[0]
            sfa=attrs["scale_factor"]
            scale_factor = sfa[0]        
            ua=attrs["units"]
            units = ua[0]
            fattrs = hdf.attributes(full=1)
            ga = fattrs["StructMetadata.0"]
            gridmeta = ga[0]

        # Construct the grid.  The needed information is in a global attribute
        # called 'StructMetadata.0'.  Use regular expressions to tease out the
        # extents of the grid.  In addition, the grid is in packed decimal
        # degrees, so we need to normalize to degrees.

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
    

    # In basemap, the sinusoidal projection is global, so we won't use it.
    # Instead we'll convert the grid back to lat/lons so we can use a local 
    # projection.
    sinu = pyproj.Proj("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
    wgs84 = pyproj.Proj("+init=EPSG:4326") 
    lon, lat= pyproj.transform(sinu, wgs84, xv, yv)

    # Apply the attributes to the data.
    invalid = np.logical_or(data < valid_range[0], data > valid_range[1])
    invalid = np.logical_or(invalid, data ==_FillValue)
    data[invalid] = np.nan
    data = (data - add_offset) /  scale_factor
    data = np.ma.masked_array(data, np.isnan(data))
    
    # A plain geographic projection looks a little warped at this scale and
    # latitude, so use a Lambert Azimuthal Equal Area projection instead.
    m = Basemap(projection='laea', resolution='l',
                lat_ts=35, lat_0=35, lon_0=-92.5,
                width=2500000, height=2000000)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(30, 45, 5), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-105, -75, 5), labels=[0, 0, 0, 1])
    m.pcolormesh(lon[::2,::2], lat[::2,::2], data[::2,::2], latlon=True)

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
    hdffile = 'MYD13A1.A2006321.h10v05.004.2006341182856.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    


