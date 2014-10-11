"""
This example code illustrates how to access and visualize a MOPITT HDF-EOS2 
swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOP02_20000303_L2V5_7_1.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

def run(FILE_NAME):
    
    with h5py.File(FILE_NAME, mode='r') as f:

        name = '/HDFEOS/SWATHS/O3NadirSwath/Data Fields/O3'
        data = f[name][:, 5]
        units = f[name].attrs['Units'].decode()
        longname = f[name].attrs['Title'].decode()
        fillvalue = f[name].attrs['_FillValue']

        data[data == fillvalue] = np.nan

        # Get the geolocation data
        latitude = f['/HDFEOS/SWATHS/O3NadirSwath/Geolocation Fields/Latitude'][:]
        longitude = f['/HDFEOS/SWATHS/O3NadirSwath/Geolocation Fields/Longitude'][:]

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 45))
    m.drawmeridians(np.arange(-180, 180, 45), labels=[True,False,False,True])
    m.scatter(longitude, latitude, c=data, s=1, cmap=plt.cm.jet,
            edgecolors=None, linewidth=0)
    m.colorbar()
    plt.title('{0} ({1})\n'.format(name, 'hPa'))
    
    fig = plt.gcf()
    plt.show()
    
    # Make an output filename out of the filename and the variable itself.
    base = os.path.splitext(os.path.basename(FILE_NAME))[0]
    varname = os.path.basename(name)
    pngfile = "{0}.{1}.png".format(os.path.basename(__file__), varname)
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'TES-Aura_L2-O3-Nadir_r0000011015_F05_07.he5'

    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
