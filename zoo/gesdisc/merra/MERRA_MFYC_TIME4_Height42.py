"""
This example code illustrates how to access and visualize a GESDISC MERRA file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MERRA_MFYC_TIME4_Height42.py

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

    DATAFIELD_NAME = 'MFYC'
    
    dset = Dataset(FILE_NAME)
    data = dset.variables[DATAFIELD_NAME][4, 42, :, :].astype(np.float64)
    
    # Replace the missing values with NaN.
    missing_value = dset.variables[DATAFIELD_NAME].missing_value
    data[data == missing_value] = np.nan
    datam = np.ma.masked_array(data, np.isnan(data))
    
    # Retrieve the geolocation data.
    latitude = dset.variables['YDim'][:]
    longitude = dset.variables['XDim'][:]
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 90., 30.))
    m.drawmeridians(np.arange(-180, 180., 45.))
    
    # Render the image in the projected coordinate system.
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, datam)
    m.colorbar()
    fig = plt.gcf()
    
    plt.title('{0} ({1})\nat TIME=4 and Height=42m'.format(
        dset.variables[DATAFIELD_NAME].long_name,
        dset.variables[DATAFIELD_NAME].units))
    plt.show()
    
    png = "{0}.{1}.png".format(os.path.basename(FILE_NAME)[:-4],
                               os.path.basename(DATAFIELD_NAME))
    fig.savefig(png)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MERRA300.prod.assim.tavg3_3d_chm_Nv.20021201.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
