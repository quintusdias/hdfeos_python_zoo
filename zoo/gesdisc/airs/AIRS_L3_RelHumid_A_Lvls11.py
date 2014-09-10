"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC AIRS grid
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AIRS_L3_RelHumid_A_Lvls11.py

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

USE_NETCDF4 = False

def run(FILE_NAME):

    # Identify the HDF-EOS2 grid data file.
    DATAFIELD_NAME = 'RelHumid_A'

    if USE_NETCDF4:
        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)

        # The variable has a fill value, 
        # so netCDF4 converts it to a float64 masked array for us.
        data = nc.variables[DATAFIELD_NAME][11,:,:]
        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]
    
    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # List available SDS datasets.
        # print hdf.datasets()

        data3D = hdf.select(DATAFIELD_NAME)
        data = data3D[11,:,:]
        # Read geolocation dataset.
        lat = hdf.select('Latitude')
        latitude = lat[:,:]
        lon = hdf.select('Longitude')
        longitude = lon[:,:]
        # Handle fill value.
        attrs = data3D.attributes(full=1)
        fillvalue=attrs["_FillValue"]
        # fillvalue[0] is the attribute value.
        fv = fillvalue[0]
        data[data == fv] = np.nan
        data = np.ma.masked_array(data, np.isnan(data))
        
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1])
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, data)
    cb = m.colorbar()
    cb.set_label('Unit:%')
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n {1} at H20PrsLvls=11'.format(basename, DATAFIELD_NAME))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.{1}.py.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'AIRS.2002.08.01.L3.RetStd_H031.v4.0.21.0.G06104133732.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
