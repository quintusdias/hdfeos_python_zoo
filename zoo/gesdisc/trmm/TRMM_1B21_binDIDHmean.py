"""
This example code illustrates how to access and visualize a GESDISC TRMM file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python TRMM_1B21_CSI_binDIDHmean_zoom.py

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

    DATAFIELD_NAME = 'binDIDHmean'
    
    dset = Dataset(FILE_NAME)
    data = dset.variables[DATAFIELD_NAME][:].astype(np.float64)
    
    # Retrieve the geolocation data.
    latitude = dset.variables['geolocation'][:,:,0]
    longitude = dset.variables['geolocation'][:,:,1]
    
    # The swath crosses the international dateline between row 8000 and 9000.
    # This causes the mesh to smear, so we'll adjust the longitude (modulus
    # 360 degrees, of course) both in the longitude array and in the basemap
    # definition to avoid that.
    longitude[longitude < -90] += 360
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-90, urcrnrlon = 270)
    
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1])
    
    # Render the image in the projected coordinate system.
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, data)
    m.colorbar()
    fig = plt.gcf()
    
    plt.title('{0}'.format(DATAFIELD_NAME))
    plt.show()
    
    pngfile = "{0}.{1}.png".format(os.path.basename(FILE_NAME[:-4]),
                                   DATAFIELD_NAME)
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = '1B21.071022.56609.6.HDF'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
