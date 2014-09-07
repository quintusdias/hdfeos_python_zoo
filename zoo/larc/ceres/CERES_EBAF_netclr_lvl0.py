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

    python CERES_EBAF_netclr_lvl0.py

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

    # Identify the data field.
    DATAFIELD_NAME = 'netclr'
    
    # Subset the data to match the size of the swath geolocation fields.
    # Turn off autoscaling, we'll handle that ourselves due to presence of
    # a valid range.
    var = nc.variables[DATAFIELD_NAME]
    data = var[0,:,:].astype(np.float64)
    latitude = nc.variables['lat'][:]
    longitude = nc.variables['lon'][:]
    
    # Apply the attributes.
    invalid = np.logical_or(data < var.valid_range[0],
                            data > var.valid_range[1])
    data[invalid] = np.nan
    datam = np.ma.masked_array(data, mask=np.isnan(data))
    
    # The data is global, so render in a global projection.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90.,90,45))
    m.drawmeridians(np.arange(-180.,180,45), labels=[True,False,False,True])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    m.colorbar()
    plt.title('{0} ({1})\n'.format(DATAFIELD_NAME, var.units))
    
    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = basename + ".png"
    fig.savefig(pngfile)
    
if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    ncfile = 'CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc.hdf'
    try:
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'], ncfile)
    except KeyError:
        fname = hdffile

    run(fname)
