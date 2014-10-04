"""
This example code illustrates how to access and visualize a MOPITT ASDC MOP03 
version 6 HDF-EOS5 grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOP02J_20131129_L2V16_2_3.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os
import re

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

def run(FILE_NAME):
    
    with h5py.File(FILE_NAME, mode='r') as f:

        group = f['/HDFEOS/GRIDS/MOP03/Data Fields']
        dsname = 'RetrievedSurfaceTemperatureDay'
        data = group[dsname][:].T
        longname = group[dsname].attrs['long_name'].decode()
        units = group[dsname].attrs['units'].decode()
        fillvalue = group[dsname].attrs['_FillValue']

        data[data == fillvalue] = np.nan
        data = np.ma.masked_array(data, np.isnan(data))

        # We could query the string dataset
        # '/HDFEOS INFORMATION/StructMetadata.0' for the geolocation
        # information, but in this case we also have lat and lon datasets.
        y = f['/HDFEOS/GRIDS/MOP03/Data Fields/Latitude'][:]
        x = f['/HDFEOS/GRIDS/MOP03/Data Fields/Longitude'][:]
        longitude, latitude = np.meshgrid(x, y)

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 45))
    m.drawmeridians(np.arange(-180, 180, 45), labels=[True,False,False,True])
    sc = m.scatter(longitude, latitude, c=data, s=1, cmap=plt.cm.jet,
                   edgecolors=None, linewidth=0)
    m.colorbar()
    plt.title('{0} ({1})\n'.format(dsname, units))
    
    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, dsname)
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOP03T-20131129-L3V4.2.1.he5'

    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)

