"""
This example code illustrates how to access and visualize a LAADS
MYD swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MYD02HKM_A2010031_0035_005_2010031183706_EV_500_RefSB_lvl0.py

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

    DATAFIELD_NAME = 'EV_500_RefSB'
    
    nc = Dataset(FILE_NAME)

    # Just read the first level, and subset the data to match the lat/lon
    # resolution.
    data = nc.variables[DATAFIELD_NAME][0,::2,::2].astype(np.float64)
    units = nc.variables[DATAFIELD_NAME].radiance_units
    long_name = nc.variables[DATAFIELD_NAME].long_name

    # The scale and offset attributes do not have standard names in this case,
    # so we have to apply the scaling equation ourselves.  Fill value is 
    # already applied, though.
    scale = nc.variables[DATAFIELD_NAME].radiance_scales[0]
    offset = nc.variables[DATAFIELD_NAME].radiance_offsets[1]
    valid_range = nc.variables[DATAFIELD_NAME].valid_range
    invalid = np.logical_or(data < valid_range[0], data > valid_range[1])
    data[invalid] = np.nan
    data = scale * (data - offset)
    data = np.ma.masked_array(data, np.isnan(data))

    latitude = nc.variables['Latitude'][:]
    longitude = nc.variables['Longitude'][:]
    
    # Use a hemispherical projection for the southern hemisphere since the
    # swath is over Antarctica.
    m = Basemap(projection='splaea', resolution='l', 
                boundinglat=-65, lon_0=0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, -50, 10), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)
    m.colorbar()
    plt.title('{0}\n({1})'.format(long_name, units))

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MYD02HKM.A2010031.0035.005.2010031183706.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
