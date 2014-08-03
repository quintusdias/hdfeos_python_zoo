"""
This example code illustrates how to access and visualize an NSIDC AMSR swath
Swath data file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AMSR_E_L2_Ocean_V06_200206190029_D_High_res_cloud.py

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
    
    nc = Dataset(FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'High_res_cloud'
    
    data = nc.variables[DATAFIELD_NAME][:].astype(np.float64)
    latitude = nc.variables['Latitude'][:]
    longitude = nc.variables['Longitude'][:]

    # There is a wrap-around effect to deal with, as some of the swath extends
    # eastward over the international dateline.  Adjust the longitude to avoid
    # the swath being smeared.
    longitude[longitude < -170] += 360

    # Apply the fill value and scaling equation.
    data[data == -9990] = np.nan
    data = data * nc.variables[DATAFIELD_NAME].Scale
    data = np.ma.masked_array(data, np.isnan(data))
    
    # Draw a polar stereographic projection using the low resolution coastline
    # database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-170, urcrnrlon=190)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 45), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180,181,45), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    m.colorbar()
    plt.title("{0} (mm)".format(DATAFIELD_NAME))

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'AMSR_E_L2_Ocean_V06_200206190029_D.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
