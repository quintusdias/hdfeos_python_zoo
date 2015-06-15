"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a NSIDC Level-2
MODIS Swath data file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD10_L2_SnowCover_P.py

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


def run():
    
    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'MOD10_L2.A2000065.0040.005.2008235221207.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Snow_Cover'
    
    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        # Subset the data to match the size of the swath geolocation fields.
        rows = slice(5, 4060, 10)
        cols = slice(5, 2708, 10)
        data = nc.variables['Snow_Cover'][rows, cols]
        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:].astype(np.float64)

        # Read geolocation dataset from HDF-EOS2 dumper output.
        GEO_FILE_NAME = 'lat_MOD10_L2.A2000065.0040.005.2008235221207.output'
        GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], 
                                     GEO_FILE_NAME)
        lat = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        latitude = lat.reshape(data.shape)
        
        GEO_FILE_NAME = 'lon_MOD10_L2.A2000065.0040.005.2008235221207.output'
        GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], 
                                     GEO_FILE_NAME)
        lon = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        longitude = lon.reshape(data.shape)

    
    # Draw a polar stereographic projection using the low resolution coastline
    # database.
    m = Basemap(projection='npstere', resolution='l',
                boundinglat=64, lon_0 = 0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(60.,81,10.))
    m.drawmeridians(np.arange(-180.,181.,30.), labels=[True,False,False,True])
    
    # Use a discretized colormap since we have only two levels.
    cmap = mpl.colors.ListedColormap(['grey','mediumblue'])
    bounds = [0, 19.5, 39]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    m.pcolormesh(longitude, latitude, data, latlon=True, cmap=cmap, norm=norm)
    
    # Must reset the alpha level to opaque for the colorbar.
    # See http://stackoverflow.com/questions/4478725/...
    # .../partially-transparent-scatter-plot-but-with-a-solid-color-bar
    color_bar = plt.colorbar()
    color_bar.set_alpha(1)
    color_bar.set_ticks([9.75, 29.25])
    color_bar.set_ticklabels(['missing data', 'ocean'])
    color_bar.draw_all()

    basename = os.path.basename(FILE_NAME)
    long_name = 'Snow Cover'
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
    
