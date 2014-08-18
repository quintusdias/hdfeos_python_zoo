"""
This example code illustrates how to access and visualize an LP_DAAC MYD
swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD11_L2_LST.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):

    DATAFIELD_NAME = 'LST'
    
    nc = Dataset(FILE_NAME)

    # Subset the data to match the swath geolocation dimension.  And no need to
    # apply the scaling equation, the netCDF4 package does this for us since
    # the scale_factor, add_offset, and non-default _FillValue attributes are
    # in place.  We will use the units and long_name attributes, though
    data = nc.variables[DATAFIELD_NAME][::5,::5]
    units = nc.variables[DATAFIELD_NAME].units
    long_name = nc.variables[DATAFIELD_NAME].long_name
    
    # Retrieve the geolocation data.
    latitude = nc.variables['Latitude'][:]
    longitude = nc.variables['Longitude'][:]
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=37.5, urcrnrlat=62.5,
                llcrnrlon=-97.5, urcrnrlon = -57.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(40, 70, 10), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-100, 60, 10), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    m.colorbar()
    plt.title('{0} ({1})'.format(long_name, units))

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MYD11_L2.A2007093.0735.005.2007101061952.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
