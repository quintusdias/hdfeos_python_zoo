"""
This example code illustrates how to access and visualize a LAADS MYD (MODIS-
AQUA) swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MYD021KM_EV_1KM_Emmisive_level0.py

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

    DATAFIELD_NAME = 'EV_1KM_Emissive'
    
    nc = Dataset(FILE_NAME)

    # Just read the first level, Band 20.
    data = nc.variables[DATAFIELD_NAME][0,:,:].astype(np.float64)
    units = nc.variables[DATAFIELD_NAME].radiance_units
    long_name = nc.variables[DATAFIELD_NAME].long_name

    # The scale and offset attributes do not have standard names in this case,
    # so we have to apply the scaling equation ourselves.
    scale = nc.variables[DATAFIELD_NAME].radiance_scales[0]
    offset = nc.variables[DATAFIELD_NAME].radiance_offsets[1]
    valid_range = nc.variables[DATAFIELD_NAME].valid_range
    fill_value = nc.variables[DATAFIELD_NAME]._FillValue
    invalid = np.logical_or(data < valid_range[0], data > valid_range[1])
    invalid = np.logical_or(invalid, data == fill_value)
    data[invalid] = np.nan
    data = scale * (data - offset)
    data = np.ma.masked_array(data, np.isnan(data))

    latitude = nc.variables['Latitude'][:]
    longitude = nc.variables['Longitude'][:]
    
    # The data is close to the equator in Africa, so a global projection is
    # not needed.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-5, urcrnrlat=30, llcrnrlon=5, urcrnrlon=45)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(0, 50, 10), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(0, 50., 10), labels=[0, 0, 0, 1])
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
    hdffile = 'MYD021KM.A2002226.0000.005.2009193222735.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
