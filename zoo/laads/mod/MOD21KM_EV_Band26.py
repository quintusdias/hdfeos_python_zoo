"""
This example code illustrates how to access and visualize a LAADS MODIS swath
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD21KM_EV_Band26.py

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

    DATAFIELD_NAME = 'EV_Band26'
    
    nc = Dataset(FILE_NAME)
    var = nc.variables[DATAFIELD_NAME]

    # Have to be very careful of the scaling equation here.
    # We'll turn autoscaling off in order to correctly scale the data.
    # Also need to subset the data to match the lat/lon dimensions.
    var.set_auto_maskandscale(False)
    data = var[2:2030:5, 2:1354:5].astype(np.double)
    data[data == var._FillValue] = np.nan
    data[data > var.valid_range[1]] = np.nan
    data = (data - var.radiance_offsets) * var.radiance_scales 
    datam = np.ma.masked_array(data, np.isnan(data))
    
    # Retrieve the geolocation data.
    longitude = nc.variables['Longitude'][:]
    latitude = nc.variables['Latitude'][:]
    
    # Render the plot in a lambert equal area projection.
    m = Basemap(projection='laea', resolution='l', lat_ts=65,
                lat_0=65, lon_0=-35,
                width=3000000,height=2500000)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(50., 91., 10.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181., 30), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    m.colorbar()
    titlestr = 'Radiance derived from Earth View Band 26 Scaled Integers\n'
    titlestr += '(watts/m^2/micrometer/steradian)'
    plt.title(titlestr)

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOD021KM.A2000055.0000.005.2010041143816.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
