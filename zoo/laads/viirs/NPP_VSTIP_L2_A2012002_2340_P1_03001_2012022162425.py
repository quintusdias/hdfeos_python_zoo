"""
This example code illustrates how to access and visualize a LAADS NPP VIIRS
swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python NPP_VSTIP_L2_A2012002_2340_P1_03001_2012022162425.py

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

    DATAFIELD_NAME = 'SurfaceTemperature'
    nc = Dataset(FILE_NAME)

    # The dataset is (6144 x 6400).  Subset it to be around than 1K x 1K
    rows = slice(0, 6144, 6)
    cols = slice(0, 6400, 6)
    data = nc.variables[DATAFIELD_NAME][rows, cols]
    
    # Retrieve the geolocation data.
    latitude = nc.variables['Latitude'][rows, cols]
    longitude = nc.variables['Longitude'][rows, cols]

    # Apply the fill value.  The valid minimum is zero, although there's no
    # attribute.
    data[data < 0] = np.nan
    data = np.ma.masked_array(data, np.isnan(data))
    
    # Render the data in a lambert azimuthal equal area projection.
    m = Basemap(projection='nplaea', resolution='l',
                boundinglat=60, lon_0=43)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(50, 90, 10), labels=[1, 0, 0, 1])
    m.drawmeridians(np.arange(-180, 180, 30))
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, data)
    m.colorbar()
    plt.title('{0}'.format(DATAFIELD_NAME))

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'NPP_VSTIP_L2.A2012002.2340.P1_03001.2012022162425.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
