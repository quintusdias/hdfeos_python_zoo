"""
This example code illustrates how to access and visualize a LaRC CERES file in
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python CER_SYN_Aqua_OTF_LTC_Sky_lvl2_Ham.py

The netCDF file must either be in your current working directory
or in a directory specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):

    nc = Dataset(FILE_NAME)

    # Subset the data for synoptic hours = 2.  
    DATAFIELD_NAME = 'LW TOA Clear-Sky'
    var = nc.variables[DATAFIELD_NAME]
    data = var[2,:,:].astype(np.float64)
    
    # The normal grid information is not present.  We have to generate the geo-
    # location data, see [1] for details.
    ysize, xsize = data.shape
    xinc = 360.0 / xsize
    yinc = 180.0 / ysize
    x0, x1 = (-180, 180)
    y0, y1 = (-90, 90)
    lon = np.linspace(x0 + xinc/2, x1 - xinc/2, xsize)
    lat = np.linspace(y0 + yinc/2, y1 - yinc/2, ysize)

    # Flip the latitude to run from 90 to -90.
    lat = lat[::-1]
    longitude, latitude = np.meshgrid(lon, lat)
    
    # The data is global, so render in a global projection.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m = Basemap(projection='hammer', lon_0=0, resolution='l')
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90.,90,45))
    m.drawmeridians(np.arange(-180.,180,45))
    m.pcolormesh(longitude, latitude, data, latlon=True)
    m.colorbar()
    plt.title('{0} at Synoptic Hours = 2\n{1}'.format(DATAFIELD_NAME, var.units))
    
    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = basename + ".png"
    fig.savefig(pngfile)
    
if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'CER_SYN_Aqua-FM3-MODIS_Edition2B_007005.20051031.hdf'
    try:
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        fname = hdffile

    run(fname)

