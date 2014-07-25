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

    python MODARNSS_EV_1KM_Emissive_level0.py

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
    
    dset = Dataset(FILE_NAME)
    var = dset.variables[DATAFIELD_NAME]

    # Have to be very careful of the scaling equation here.
    # We'll turn autoscaling off in order to correctly scale the data.
    var.set_auto_maskandscale(False)
    data = var[0,:,:].astype(np.double)
    data[data == var._FillValue] = np.nan
    data[data > var.valid_range[1]] = np.nan
    data = (data - var.radiance_offsets[0]) * var.radiance_scales[0] 
    datam = np.ma.masked_array(data, np.isnan(data))
    
    # Retrieve the geolocation data.
    longitude = dset.variables['Longitude'][:]
    latitude = dset.variables['Latitude'][:]
    
    # Render the plot in a cylindrical projection.
    m = Basemap(projection='cyl', resolution='l', 
                llcrnrlat=-12, urcrnrlat = -9,
                llcrnrlon=-64, urcrnrlon = -61)

    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-12., -8., 1.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-64, -60., 1), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    cb = m.colorbar()
    titlestr = '{0}'.format(var.long_name)
    plt.title(titlestr)

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MODARNSS.Abracos_Hill.A2000080.1515.005.2007164153544.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
