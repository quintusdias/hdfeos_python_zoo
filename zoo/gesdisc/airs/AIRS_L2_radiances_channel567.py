"""
This example code illustrates how to access and visualize a GESDISC AIRS grid
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AIRS_L2_radiances_channel567.py

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

USE_NETCDF4 = False

def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'AIRS.2002.12.31.001.L2.CC_H.v4.0.21.0.G06100185050.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the HDF-EOS2 swath data file.
    DATAFIELD_NAME = 'radiances'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        data = nc.variables['radiances'][:, :, 567].astype(np.float64)

        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data3D = hdf.select(DATAFIELD_NAME)
        data = data3D[:, :, 567]

        # Read geolocation dataset.
        lat = hdf.select('Latitude')
        latitude = lat[:]
        lon = hdf.select('Longitude')
        longitude = lon[:]

    # Replace the filled value with NaN, replace with a masked array.
    data[data == -9999] = np.nan
    datam = np.ma.masked_array(data, np.isnan(data))

    # Draw a polar stereographic projection using the low resolution coastline
    # database.
    m = Basemap(projection='spstere', resolution='l',
                boundinglat=-65, lon_0=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-80., -50., 5.))
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[1, 0, 0, 1])

    # Render the image in the projected coordinate system.
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, datam)
    m.colorbar()

    # See page 101 of "AIRS Version 5.0 Released Files Description"
    # document [1] for unit specification.
    units = 'mW/m**2/cm**-1/sr'
    plt.title('{0} ({1}) at channel 567'.format(DATAFIELD_NAME, units))

    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig = plt.gcf()
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
