"""
This example code illustrates how to access and visualize a NSIDC MOD29
Level 2 HDF-EOS2 Swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD29_A2013196_1250_005_2013196195940_hdf.py

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
    DATAFIELD_NAME = 'Ice_Surface_Temperature'
    
    # Subset the data to match the size of the swath geolocation fields.
    rows = slice(2, 2030, 5)
    cols = slice(2, 1354, 5)
    data = nc.variables[DATAFIELD_NAME][rows, cols]
    latitude = nc.variables['Latitude'][:]
    longitude = nc.variables['Longitude'][:]
    
    # Retrieve any attributes we might need.
    units = nc.variables[DATAFIELD_NAME].units
    scale = nc.variables[DATAFIELD_NAME].scale_factor
    fillvalue = nc.variables[DATAFIELD_NAME]._FillValue
    valid_range = nc.variables[DATAFIELD_NAME].valid_range
    
    # The scale factor and add_offset are already applied by netCDF4 when the
    # data was read.  We still have to restrict the data to the valid range,
    # though.  Since the array will have NaNs, we then need to turn it into
    # a masked array in order to render it.
    data[data < float(valid_range[0])*scale] = np.nan
    data[data > float(valid_range[1])*scale] = np.nan
    datam = np.ma.masked_array(data, mask=np.isnan(data))
    
    # Draw a southern polar stereographic projection using the low resolution
    # coastline database.
    m = Basemap(projection='spstere', resolution='l',
                boundinglat=-64, lon_0 = 0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-80.,-59,10.))
    m.drawmeridians(np.arange(-180.,179.,30.), labels=[True,False,False,True])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    m.colorbar()
    plt.title('{0} ({1})\n'.format(DATAFIELD_NAME, units))
    
    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = basename + ".png"
    fig.savefig(pngfile)
    
if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOD29.A2013196.1250.005.2013196195940.hdf'
    try:
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        fname = hdffile

    run(fname)
