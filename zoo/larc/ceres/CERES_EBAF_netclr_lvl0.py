"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LaRC CERES file in
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python CERES_EBAF_netclr_lvl0.py

The netCDF file must either be in your current working directory
or in a directory specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = False

def run(FILE_NAME):

    # Identify the data field.
    DATAFIELD_NAME = 'netclr'

    if USE_NETCDF4:
        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)
    
        # Subset the data to match the size of the swath geolocation fields.
        # Turn off autoscaling, we'll handle that ourselves due to presence of
        # a valid range.
        var = nc.variables[DATAFIELD_NAME]
        data = var[0,:,:].astype(np.float64)
        latitude = nc.variables['lat'][:]
        longitude = nc.variables['lon'][:]

        # Read attributes.
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')]
        units = ncvar.units
        long_name = ncvar.long_name

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)
        
        # Read dataset.
        data3D = hdf.select(DATAFIELD_NAME)
        data = data3D[0,:,:]

        # Read geolocation datasets.
        lat = hdf.select('lat')
        latitude = lat[:]
        lon = hdf.select('lon')
        longitude = lon[:]

        # Read attributes.
        attrs = data3D.attributes(full=1)
        ua=attrs["units"]
        units = ua[0]
        vra=attrs["valid_range"]
        valid_range = vra[0]
        lna=attrs["long_name"]
        long_name = lna[0]    

    # Apply the valid_range attribute.
    invalid = np.logical_or(data < valid_range[0],
                            data > valid_range[1])
    data[invalid] = np.nan
    datam = np.ma.masked_array(data, mask=np.isnan(data))
    
    # The data is global, so render in a global projection.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1} at time=0'.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)
    
if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    ncfile = 'CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc.hdf'
    try:
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'], ncfile)
    except KeyError:
        fname = hdffile

    run(fname)
