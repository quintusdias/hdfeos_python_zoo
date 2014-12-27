"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an NSIDC AMSR_E 
version 3 L2A HDF-EOS2 swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AMSR_E_L2A_BrightnessTemperatures_V12_201110032238_D_hdf.py

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

    DATAFIELD_NAME = '89.0V_Res.5B_TB_(not-resampled)'

    if USE_NETCDF4:
        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)

        data = nc.variables[DATAFIELD_NAME][:].astype(np.float64)
        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]
    
        # Replace the filled value with NaN, replace with a masked array.
        # Apply the scaling equation.  These attributes are named in a VERY
        # non-standard manner.
        scale_factor = getattr(nc.variables[DATAFIELD_NAME], 'SCALE FACTOR')
        add_offset = nc.variables[DATAFIELD_NAME].OFFSET

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:,:].astype(np.float64)

        # Read geolocation dataset.
		# This product has multiple 'Latitude' and 'Longitude' pair under different groups.
        lat = hdf.select(hdf.reftoindex(192)) # Use HDFView to get ref number 192.
        latitude = lat[:,:]
        lon = hdf.select(hdf.reftoindex(194)) # Use HDFView to get ref number 194.
        longitude = lon[:,:]

        # Retrieve attributes.
        attrs = data2D.attributes(full=1)
        sfa=attrs["SCALE FACTOR"]
        scale_factor = sfa[0]        
        aoa=attrs["OFFSET"]
        add_offset = aoa[0]
        
    data[data == -32768] = np.nan
    data = data * scale_factor + add_offset
    datam = np.ma.masked_array(data, np.isnan(data))

    units = "degrees K"
    long_name = DATAFIELD_NAME

    # Since the swath starts near the south pole, but also extends over the
    # north pole, the equidistant cylindrical becomes a possibly poor choice
    # for a projection.  We show the full global map plus a limited polar map.
    fig = plt.figure(figsize=(15, 6))
    ax1 = plt.subplot(1, 2, 1)
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181., 45), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)

    ax2 = plt.subplot(1, 2, 2)
    m = Basemap(projection='npstere', resolution='l',
                boundinglat=65, lon_0=0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(60, 81, 10), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 30.), labels=[1, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)

    cax = plt.axes([0.92, 0.1, 0.03, 0.8])
    cb = plt.colorbar(cax=cax)
    cb.set_label(units)
    

    basename = os.path.basename(FILE_NAME)
    fig = plt.gcf()
    fig.suptitle('{0}\n{1}'.format(basename, long_name))
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'AMSR_E_L2A_BrightnessTemperatures_V12_201110032238_D.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
