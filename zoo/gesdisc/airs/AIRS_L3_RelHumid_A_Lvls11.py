"""
This example code illustrates how to access and visualize a GESDISC_AIRS
Grid in Python.  If you have any questions, suggestions, comments  on
this example, please use the HDF-EOS Forum (http://hdfeos.org/forums).
If you would like to see an  example of any other NASA HDF/HDF-EOS data
product that is not listed in the HDF-EOS Comprehensive Examples page
(http://hdfeos.org/zoo), feel free to contact us at eoshelp@hdfgroup.org
or post it at the HDF-EOS Forum (http://hdfeos.org/forums).
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):

    # Identify the HDF-EOS2 grid data file.
    DATAFIELD_NAME = 'RelHumid_A'
    
    dset = Dataset(FILE_NAME)

    # The variable has a fill value, so netCDF4 converts it to a float64 masked
    # array for us.
    data = dset.variables[DATAFIELD_NAME][11,:,:]
    
    latitude = dset.variables['Latitude'][:]
    longitude = dset.variables['Longitude'][:]
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 90., 30.))
    m.drawmeridians(np.arange(-180., 181., 45.))
    
    # Render the image in the projected coordinate system.
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, data)
    m.colorbar()
    fig = plt.gcf()
    
    plt.title('{0} at H20PrsLvls=11'.format(DATAFIELD_NAME, ""))
    plt.show()
    
    pngfile = "{0}.{1}.png".format(os.path.basename(FILE_NAME[:-4]),
                                   DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'AIRS.2002.08.01.L3.RetStd_H031.v4.0.21.0.G06104133732.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
