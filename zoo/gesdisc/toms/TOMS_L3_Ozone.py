"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC TOMS grid
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python TOMS_L3_Ozone.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_PYHDFEOS = True
USE_NETCDF4 = False
USE_PYHDF = False

FILENAME = 'TOMS-EP_L3-TOMSEPL3_2000m0101_v8.HDF'

def run():

    DATAFIELD_NAME = 'Ozone'

    path = fullpath(FILENAME)

    if USE_PYHDFEOS:

        from pyhdfeos import GridFile    
        grid = GridFile(path).grids['TOMS Level 3']
        data = grid.fields[DATAFIELD_NAME][:]
        missing_value = grid.fields[DATAFIELD_NAME].attrs['missing_value']
        units = grid.fields[DATAFIELD_NAME].attrs['units']
        long_name = grid.fields[DATAFIELD_NAME].attrs['long_name']

        latitude, longitude = grid[:]
        #latitude = nc.variables['YDim:TOMS Level 3'][:]
        #longitude = nc.variables['XDim:TOMS Level 3'][:]

    elif USE_NETCDF4:

        from netCDF4 import Dataset    
        nc = Dataset(path)
        data = nc.variables[DATAFIELD_NAME][:]
        missing_value = nc.variables[DATAFIELD_NAME].missing_value
        units = nc.variables[DATAFIELD_NAME].units
        long_name = nc.variables[DATAFIELD_NAME].long_name

        latitude = nc.variables['YDim:TOMS Level 3'][:]
        longitude = nc.variables['XDim:TOMS Level 3'][:]

    elif USE_PYHDF:

        from pyhdf.SD import SD, SDC
        hdf = SD(path, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:,:]

        # Retrieve the attributes.
        attrs = data2D.attributes(full=1)
        mva=attrs["missing_value"]
        missing_value = mva[0]
        lna=attrs["long_name"]
        long_name = lna[0]
        ua=attrs["units"]
        units = ua[0]        

        # Read geolocation dataset.
        lat = hdf.select('YDim:TOMS Level 3')
        latitude = lat[:]
        lon = hdf.select('XDim:TOMS Level 3')
        longitude = lon[:]

    # Replace the missing values with NaN.        
    data[data == missing_value] = np.nan
    datam = np.ma.masked_array(data, np.isnan(data))

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(path)
    plt.title('{0}\n{1}'.format(basename, long_name))

    fig = plt.gcf()
    # plt.show()
    
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


def fullpath(filename):
    try:
        # assume it's NOT in the current directory, but rather at a place
        # pointed to by this particular environment variable.
        path = os.path.join(os.environ['HDFEOS_ZOO_DIR'], filename)
    except KeyError:
        # assume it's in the current directory
        path = filename
    return path

if __name__ == "__main__":
    run()
    
