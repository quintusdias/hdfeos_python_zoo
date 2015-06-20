"""
This example code illustrates how to access and visualize a PO.DAAC AVHRR
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python PODAAC_L3_bsst.py

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

USE_NETCDF4 = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = '2006001-2006005.s0454pfrt-bsst.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'bsst'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)

        # Read the dataset and geolocation data.
        # Subset by a factor of 8 to speed up plotting.
        # Turn off autoscaling, we'll handle that ourselves due to non-standard
        # naming of the offset attribute.
        var = nc.variables[DATAFIELD_NAME]
        var.set_auto_maskandscale(False)
        data = var[::8, ::8].astype(np.float64)
        latitude = nc.variables['lat'][::8]
        longitude = nc.variables['lon'][::8]

        scale_factor = var.scale_factor
        add_offset = var.add_off
        units = var.units

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read the dataset and geolocation data.
        # Subset by a factor of 8 to speed up plotting.
        variable = hdf.select(DATAFIELD_NAME)
        data = variable[::8, ::8].astype(np.float64)
        latitude = hdf.select('lat')[::8]
        longitude = hdf.select('lon')[::8]

        # Retrieve attributes
        attrs = variable.attributes(full=1)
        scale_factor = attrs["scale_factor"][0]
        add_offset = attrs["add_off"][0]
        units = attrs["units"][0]

    # Apply the attributes.  By inspection, fill value is 0
    data[data == 0] = np.nan
    data = data * scale_factor + add_offset
    datam = np.ma.masked_array(data, mask=np.isnan(data))

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
