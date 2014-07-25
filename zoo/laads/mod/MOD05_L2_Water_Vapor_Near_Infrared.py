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

    python MOD05_L2_Water_Vapor_Near_Infrared.py

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

    DATAFIELD_NAME = 'Water_Vapor_Near_Infrared'
    
    dset = Dataset(FILE_NAME)
    var = dset.variables[DATAFIELD_NAME]

    # The scaling equation to be used here is not 
    #
    #     data = data * scale + offset
    #
    # We'll turn autoscaling off in order to correctly scale the data.
    var.set_auto_maskandscale(False)
    data = var[:].astype(np.double)
    data[data < var.valid_range[0]] = np.nan
    data[data > var.valid_range[1]] = np.nan
    data[data == var._FillValue] = np.nan
    data = (data - var.add_offset) * var.scale_factor 
    datam = np.ma.masked_array(data, np.isnan(data))
    
    # Retrieve the geolocation data.
    longitude = dset.variables['Longitude'][:]
    latitude = dset.variables['Latitude'][:]
    
    # Render the plot in a south plar stereographic projection.
    m = Basemap(projection='spstere', resolution='l',
                boundinglat=-60, lon_0=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 50., 10.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181., 30), labels=[0, 0, 0, 1])
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, datam)
    cb = m.colorbar()
    plt.title(var.long_name)

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOD05_L2.A2010001.0000.005.2010005211557.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
