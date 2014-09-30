"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LAADS MODIS swath
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MODATML2_Cloud_Fraction.py

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

    DATAFIELD_NAME = 'Cloud_Fraction'
    if USE_NETCDF4:        
        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)
        var = nc.variables[DATAFIELD_NAME]

        # Have to be very careful of the scaling equation here.
        # We'll turn autoscaling off in order to correctly scale the data.
        var.set_auto_maskandscale(False)
        data = var[:].astype(np.double)

        # Retrieve the geolocation data.
        longitude = nc.variables['Longitude'][:]
        latitude = nc.variables['Latitude'][:]

        # Retrieve attributes.
        scale_factor = var.scale_factor 
        add_offset = var.add_offset
        _FillValue = var._FillValue
        valid_min = var.valid_range[0]
        valid_max = var.valid_range[1]
        long_name = var.long_name
        units = var.units

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:,:].astype(np.double)

        # Read geolocation dataset.
        lat = hdf.select('Latitude')
        latitude = lat[:,:]
        lon = hdf.select('Longitude')
        longitude = lon[:,:]
        
        # Retrieve attributes.
        attrs = data2D.attributes(full=1)
        lna=attrs["long_name"]
        long_name = lna[0]
        aoa=attrs["add_offset"]
        add_offset = aoa[0]
        fva=attrs["_FillValue"]
        _FillValue = fva[0]
        sfa=attrs["scale_factor"]
        scale_factor = sfa[0]        
        vra=attrs["valid_range"]
        valid_min = vra[0][0]        
        valid_max = vra[0][1]        
        ua=attrs["units"]
        units = ua[0]

        # Retrieve attributes for lat/lon.
        attrs = lat.attributes(full=1)
        aoa=attrs["add_offset"]
        lat_add_offset = aoa[0]
        sfa=attrs["scale_factor"]
        lat_scale_factor = sfa[0]        

        attrs = lon.attributes(full=1)
        aoa=attrs["add_offset"]
        lon_add_offset = aoa[0]
        sfa=attrs["scale_factor"]
        lon_scale_factor = sfa[0]        

        # NetCDF does the following automatically but PyHDF doesn't
        latitude = (latitude - lat_add_offset) * lat_scale_factor
        longitude = (longitude - lon_add_offset) * lon_scale_factor


    invalid = np.logical_or(data > valid_max,
                            data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = (data - add_offset) * scale_factor 
    data = np.ma.masked_array(data, np.isnan(data))
    
    # Render the plot in a lambert equal area projection.
    m = Basemap(projection='laea', resolution='l', lat_ts=63,
                lat_0=63, lon_0=-45,
                width=1500000,height=1000000)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(50., 90., 10), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-55, -25., 10), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)

    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}\n'.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MODATML2.A2000055.0000.005.2006253045900.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
