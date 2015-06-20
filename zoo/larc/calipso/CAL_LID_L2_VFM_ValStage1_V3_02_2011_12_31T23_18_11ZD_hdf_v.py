"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LaRC CALIPSO file
in file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python CAL_LID_L2_VFM-ValStage1-V3-02.2011-12-31T23-18-11ZD.hdf.v.py

The netCDF file must either be in your current working directory
or in a directory specified by the environment variable HDFEOS_ZOO_DIR.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import colors

USE_NETCDF4 = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'CAL_LID_L2_VFM-ValStage1-V3-02.2011-12-31T23-18-11ZD.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Feature_Classification_Flags'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)

        # Subset the data to match the size of the swath geolocation fields.
        # Turn off autoscaling, we'll handle that ourselves due to presence of
        # a valid range.
        data = nc.variables[DATAFIELD_NAME][:]

        # Read geolocation datasets.
        lat = nc.variables['Latitude'][:]
        lat = np.squeeze(lat)

    else:

        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:]

        # Read geolocation datasets.
        latitude = hdf.select('Latitude')
        lat = latitude[:]

    # Extract Feature Type only (1-3 bits) through bitmask.
    data = data & 7

    # Subset latitude values for the region of interest (40N to 62N).
    # See the output of the
    # CAL_LID_L2_VFM-ValStage1-V3-02.2011-12-31T23-18-11ZD.hdf.py example.
    lat = lat[3500:4000]
    size = lat.shape[0]

    # You can visualize other blocks by changing subset parameters.
    # data2d = data[3500:3999, 0:164]    # 20.2km to 30.1km
    # data2d = data[3500:3999, 165:1164] #  8.2km to 20.2km

    # data2d = data[3500:4000, 1165:]  # -0.5km to  8.2km
    data2d = data[3500:4000, 1165:]  # -0.5km to  8.2km
    data3d = np.reshape(data2d, (size, 15, 290))
    data = data3d[:, 0, :]

    # Focus on cloud (=2) data only.
    data[data > 2] = 0
    data[data < 2] = 0
    data[data == 2] = 1

    # Generate altitude data according to file specification [1].
    alt = np.zeros(290)

    # You can visualize other blocks by changing subset parameters.
    #  20.2km to 30.1km
    # for i in range (0, 54):
    #       alt[i] = 20.2 + i*0.18;

    #  8.2km to 20.2km
    # for i in range (0, 199):
    #       alt[i] = 8.2 + i*0.06;

    # -0.5km to 8.2km
    for i in range(0, 289):
        alt[i] = -0.5 + i * 0.03

    # Contour the data on a grid of longitude vs. pressure
    latitude, altitude = np.meshgrid(lat, alt)

    # Make a color map of fixed colors.
    cmap = colors.ListedColormap(['white', 'blue', 'blue'])

    # Define the bins and normalize.
    bounds = np.linspace(0, 2, 3)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    long_name = 'Feature Type (Bits 1-3) in Feature Classification Flag'
    basename = os.path.basename(FILE_NAME)
    plt.contourf(latitude, altitude, np.rot90(data, 3), cmap=cmap)
    plt.title('{0}\n{1}'.format(basename, long_name))
    plt.xlabel('Latitude (degrees north)')
    plt.ylabel('Altitude (km)')

    fig = plt.gcf()

    # Create a second axes for the discrete colorbar.
    ax2 = fig.add_axes([0.93, 0.2, 0.01, 0.6])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, boundaries=bounds)
    cb.ax.set_yticklabels(['Others', '', 'Cloud'], fontsize=8)

    plt.show(block=False)
    pngfile = "{0}.v.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
