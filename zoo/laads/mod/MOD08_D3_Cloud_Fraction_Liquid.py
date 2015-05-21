"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LAADS
MOD08 grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD08_D3_Cloud_Fraction_Liquid.py

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
USE_NETCDF = True

def run(FILE_NAME):
    
    # Identify the data field.
    GRID_NAME = 'mod08'
    DATAFIELD_NAME = 'Cloud_Fraction_Liquid'

    if USE_GDAL:
        import gdal
        gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                         GRID_NAME,
                                                         DATAFIELD_NAME)
        gdset = gdal.Open(gname)
        data = gdset.ReadAsArray().astype(np.float64)

        # Read parameters for constructing the grid.
        x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
        nx, ny = (gdset.RasterXSize, gdset.RasterYSize)


        # Read fill value, valid range, scale factor, add_offset attributes.
        metadata = gdset.GetMetadata()
        valid_range = [float(x) for x in metadata['valid_range'].split(', ')]
        _FillValue = float(metadata['_FillValue'])
        scale_factor = float(metadata['scale_factor'])                   
        add_offset = float(metadata['add_offset'])
        units = metadata['units']
        long_name = metadata['long_name']
        del gdset

        x1 = x0 + xinc * nx
        y1 = y0 + yinc * ny
                           
    elif USE_NETCDF:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        ncvar = nc.variables[DATAFIELD_NAME]
        ncvar.set_auto_maskandscale(False)
        data = ncvar[:].astype(np.float64)

        valid_range = ncvar.valid_range
        _FillValue = ncvar._FillValue
        units = ncvar.units
        long_name = ncvar.long_name
        scale_factor = ncvar.scale_factor
        add_offset = ncvar.add_offset

        # The geolocation information is in a global attribute called 
        # 'StructMetadata.0'  Use regular expressions to tease out the extents
        # of the grid.
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

        # Convert to decimal degrees.
        x0 = x0 / 1000000
        y0 = y0 / 1000000
        x1 = x1 / 1000000
        y1 = y1 / 1000000
        
        ny, nx = data.shape

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:,:].astype(np.double)
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
        
	    # This product uses geographic projection.  Ideally, these
	    # parameters should be obtained by parsing StructMetadata
	    # attribute but we assume that user has checked it with
	    # HDFView.  Upper left corner: HDF-EOS2 convention
        x0 = -180 
        y0 = 90
        # Grid spacing
        xinc = 1
        yinc = -1     
        # Grid size
        nx = 360
        ny = 180

        x1 = x0 + xinc * nx
        y1 = y0 + yinc * ny
                           
    x = np.linspace(x0, x1, nx)
    y = np.linspace(y0, y1, ny)
    lon, lat = np.meshgrid(x, y)

    import pdb; pdb.set_trace()
    invalid = data < valid_range[0]
    invalid = np.logical_or(invalid, data > valid_range[1])
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan

    data =  scale_factor * (data - add_offset)
    data = np.ma.masked_array(data, np.isnan(data))
    
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}\n'.format(basename, long_name), fontsize=11)
    
    fig = plt.gcf()
    # plt.show()
    
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)




if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOD08_D3.A2010001.005.2010006233008.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    

