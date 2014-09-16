"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans


This example code illustrates how to access and visualize a GESDISC MERRA file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MERRA_PLE_TIME1_Height72.py

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

    DATAFIELD_NAME = 'PLE'
    
    if USE_NETCDF4:
        from netCDF4 import Dataset    
        nc = Dataset(FILE_NAME)
        data = nc.variables[DATAFIELD_NAME][0, 72, :, :].astype(np.float64)
    
        # Retrieve the attributes.
        missing_value = nc.variables[DATAFIELD_NAME].missing_value
        long_name = nc.variables[DATAFIELD_NAME].long_name
        units = nc.variables[DATAFIELD_NAME].units

        # Retrieve the geolocation data.
        latitude = nc.variables['YDim'][:]
        longitude = nc.variables['XDim'][:]

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data4D = hdf.select(DATAFIELD_NAME)
        data = data4D[0,72,:,:].astype(np.float64)

        # Retrieve the attributes.
        attrs = data4D.attributes(full=1)
        mva=attrs["missing_value"]
        missing_value = mva[0]
        lna=attrs["long_name"]
        long_name = lna[0]
        ua=attrs["units"]
        units = ua[0]        

        # Read geolocation dataset.
        lat = hdf.select('YDim')
        latitude = lat[:]
        lon = hdf.select('XDim')
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
    cb.set_label(units,  labelpad=-40, y=1.05)


    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1} at TIME=1 and Height=72'.format(basename,long_name))

    fig = plt.gcf()
    # plt.show()

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MERRA300.prod.assim.inst3_3d_chm_Ne.20021201.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
