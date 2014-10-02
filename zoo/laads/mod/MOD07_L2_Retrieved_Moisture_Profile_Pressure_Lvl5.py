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

    python MOD07_L2_Retrieved_Moisture_Profile_Pressure_Lvl5.py

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

    DATAFIELD_NAME = 'Retrieved_Moisture_Profile'

    if USE_NETCDF4:        
        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)
        var = nc.variables[DATAFIELD_NAME]

        # The scaling equation to be used here is not 
        #
        #     data = data * scale + offset
        #
        # We'll turn autoscaling off in order to correctly scale the data.
        var.set_auto_maskandscale(False)
        data = var[5,:,:].astype(np.double)
    
        # Retrieve the geolocation data.
        longitude = nc.variables['Longitude'][:]
        latitude = nc.variables['Latitude'][:]

        # Retrieve attributes.
        scale_factor = var.scale_factor 
        add_offset = var.add_offset
        _FillValue = var._FillValue
        long_name = var.long_name
        units = var.units

        # Retrieve dimension name.
        dimname = var.dimensions[0]

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data3D = hdf.select(DATAFIELD_NAME)
        data = data3D[5,:,:].astype(np.double)

        # Read geolocation dataset.
        lat = hdf.select('Latitude')
        latitude = lat[:,:]
        lon = hdf.select('Longitude')
        longitude = lon[:,:]

        # Retrieve attributes.
        attrs = data3D.attributes(full=1)
        lna=attrs["long_name"]
        long_name = lna[0]
        aoa=attrs["add_offset"]
        add_offset = aoa[0]
        fva=attrs["_FillValue"]
        _FillValue = fva[0]
        sfa=attrs["scale_factor"]
        scale_factor = sfa[0]        
        ua=attrs["units"]
        units = ua[0]

        # Retrieve dimension name.
        dim = data3D.dim(0)
        dimname = dim.info()[0]



    data[data == _FillValue] = np.nan
    data = (data - add_offset) * scale_factor 
    datam = np.ma.masked_array(data, np.isnan(data))
    
    # Render the plot in a south plar stereographic projection.
    m = Basemap(projection='spstere', resolution='l',
                boundinglat=-60, lon_0=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 50., 10.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181., 30), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1} at {2}=5'.format(basename, long_name, dimname))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOD07_L2.A2010001.0000.005.2010004001518.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
