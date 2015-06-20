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

    python CAL_LID_L2_VFM_ValStage1_V3_02_2011_12_31T23_18_11ZD_hdf.py

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
        data = nc.variables[DATAFIELD_NAME][:, 1256]

        # Read geolocation datasets.
        lat = nc.variables['Latitude'][:]
        lon = nc.variables['Longitude'][:]

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:, 1256]

        # Read geolocation datasets.
        latitude = hdf.select('Latitude')
        lat = latitude[:]
        longitude = hdf.select('Longitude')
        lon = longitude[:]

    # Subset data. Otherwise, all points look black.
    lat = lat[::10]
    lon = lon[::10]
    data = data[::10]

    # Extract Feature Type only through bitmask.
    data = data & 7

    # Make a color map of fixed colors.
    cmap = colors.ListedColormap(['black', 'blue', 'yellow', 'green', 'red',
                                  'purple', 'gray', 'white'])

    # The data is global, so render in a global projection.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 90, 45))
    m.drawmeridians(np.arange(-180, 180, 45),
                    labels=[True, False, False, True])
    x, y = m(lon, lat)
    i = 0
    for feature in data:
        m.plot(x[i], y[i], 'o', color=cmap(feature),  markersize=3)
        i = i + 1

    long_name = 'Feature Type at Altitude = 2500m'
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, long_name))

    fig = plt.gcf()

    # define the bins and normalize
    bounds = np.linspace(0, 8, 9)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # create a second axes for the colorbar
    ax2 = fig.add_axes([0.93, 0.2, 0.01, 0.6])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
                                   spacing='proportional', ticks=bounds,
                                   boundaries=bounds, format='%1i')

    yticklabels = ['invalid', 'clear', 'cloud', 'aerosol', 'strato', 'surface',
                   'subsurf', 'no signal']
    cb.ax.set_yticklabels(yticklabels, fontsize=5)

    plt.show(block=False)
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
