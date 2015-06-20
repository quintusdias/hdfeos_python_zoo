"""
This example code illustrates how to access and visualize a PO.DAAC QuikSCAT
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python QS_XWGRDS_des_avg_wind_speed.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4=False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'QS_XWGRD3_2008001.20080021608.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'des_avg_wind_speed'
    
    if USE_NETCDF4:

        from netCDF4 import Dataset
        nc = Dataset(FILE_NAME)
    
        # Turn off autoscaling, as the fillvalue attribute is not set.
        var = nc.variables[DATAFIELD_NAME]
        var.set_auto_maskandscale(False)
        data = var[:]
    
        # Retrieve the needed attributes.  By inspection, the fill value is 0.
        fillvalue = 0
        scale = var.scale_factor
        offset = var.add_offset
        units = var.units

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read the dataset and geolocation data.
        variable = hdf.select(DATAFIELD_NAME)
        data = variable[:].astype(np.float64)

        # Retrieve the needed attributes.  By inspection, the fill value is 0.
        fillvalue = 0
        attrs = variable.attributes(full=1)
        scale = attrs["scale_factor"][0]
        offset = attrs["add_offset"][0]
        units = attrs["units"][0]

    invalid = data == fillvalue
    data = data * scale + offset
    data[invalid] = np.nan
    datam = np.ma.masked_array(data, np.isnan(data))

    # Calculate lat and lon according to [1]
    latdim, londim = data.shape
    latinc = 180 / latdim
    loninc = 360 / londim
    x = np.linspace(0 + loninc/2, 360 - loninc/2, londim)
    y = np.linspace(-90 + latinc/2, 90 - latinc/2, latdim)
    longitude, latitude = np.meshgrid(x, y)

    # Draw a southern polar stereographic projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 45))
    m.drawmeridians(np.arange(-180, 180, 45),
                    labels=[True, False, False, True])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    m.colorbar()
    plt.title('{0} ({1})\n'.format(DATAFIELD_NAME, units))

    fig = plt.gcf()
    plt.show(block=False)

    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = basename + ".png"
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
