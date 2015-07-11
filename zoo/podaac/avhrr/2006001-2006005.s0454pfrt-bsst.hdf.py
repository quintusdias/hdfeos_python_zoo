"""
Copyright (C) 2015 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a PO.DAAC AVHRR
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python 2006001-2006005.s0454pfrt-bsst.hdf.py

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

    # Identify the data field.
    DATAFIELD_NAME = 'bsst'

    if USE_NETCDF4:
        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)
        # Subset the data to match the size of the swath geolocation fields.
        # Turn off autoscaling, we'll handle that ourselves due to non-standard
        # naming of the offset attribute.
        var = nc.variables[DATAFIELD_NAME]
        var.set_auto_maskandscale(False)
        lat = nc.variables['lat']
        lon = nc.variables['lon']
    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        var = hdf.select(DATAFIELD_NAME)
        lat = hdf.select('lat')
        lon = hdf.select('lon')

    latitude = lat[::8]
    longitude = lon[::8]
    data = var[::8, ::8].astype(np.float64)
    
    # Apply the attributes.  By inspection, fill value is 0
    data[data==0] = np.nan
    data = data * var.scale_factor + var.add_off
    datam = np.ma.masked_array(data, mask=np.isnan(data))
    
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 45))
    m.drawmeridians(np.arange(-180, 180, 45), labels=[True,False,False,True])
    m.pcolormesh(longitude, latitude, datam, latlon=True)

    cax = plt.axes([0.92, 0.3, 0.01, 0.4])
    cb = plt.colorbar(cax=cax)
    units = 'degrees-C'
    cb.set_label(units)    
    
    basename = os.path.basename(FILE_NAME)
    fig = plt.gcf()
    # plt.show()
    long_name = 'Sea Surface Temperature ('+DATAFIELD_NAME+')'
    fig.suptitle('{0}\n{1}'.format(basename, long_name))
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)
    
if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = '2006001-2006005.s0454pfrt-bsst.hdf'
    try:
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        fname = hdffile

    run(fname)
