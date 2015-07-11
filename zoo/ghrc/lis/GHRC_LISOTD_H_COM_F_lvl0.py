"""
Copyright (C) 2015 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GHRC HDF4 file in
Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python GHRC_LISOTD_H_COM_F_lvl0.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.

Last Update: 2015/06/08
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
USE_NETCDF4 = False

def run(FILE_NAME):

    # Identify the data field.
    DATAFIELD_NAME = 'HRAC_COM_FR'

    if USE_NETCDF4:
        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)
        # Subset the data to match the size of the swath geolocation fields.
        # Turn off autoscaling, we'll handle that ourselves due to the existance
        # of the valid range attribute.
        var = nc.variables[DATAFIELD_NAME]
        _FillValue = var._FillValue
        units = var.units
        long_name = var.long_name
        var.set_auto_maskandscale(False)
        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]
    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        var = hdf.select(DATAFIELD_NAME)

        # Retrieve attributes.
        attrs = var.attributes(full=1)
        fva=attrs["_FillValue"]
        _FillValue = fva[0]   
        ua=attrs["units"]
        units = ua[0]
        lna=attrs["long_name"]
        long_name = lna[0]

        lat = hdf.select('Latitude')
        lon = hdf.select('Longitude')
        latitude = lat[:]
        longitude = lon[:]
    
    data = var[:,:,0]
    data[data == _FillValue] = np.nan
    datam = np.ma.masked_array(data, np.isnan(data))


    # Retrieve the geolocation.  There's a minor wraparound issue.
    latitude[0] += 180
    longitude[0] += 360

    latitude -= 90
    longitude -= 180
    
    # Draw a southern polar stereographic projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 45))
    m.drawmeridians(np.arange(-180, 180, 45), labels=[True,False,False,True])
    m.pcolormesh(longitude, latitude, datam, latlon=True)

    cb = m.colorbar()
    cb.set_label(units)    

    basename = os.path.basename(FILE_NAME)
    long_name = long_name + ' at Day of year=0'
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)
    
if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'LISOTD_HRAC_V2.2.hdf'
    try:
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        fname = hdffile

    run(fname)
