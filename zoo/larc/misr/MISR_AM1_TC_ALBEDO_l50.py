"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an LaRC MISR SOM
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MISR_AM1_TC_ALBEDO_l50.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np
from pyhdf.SD import SD, SDC


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'MISR_AM1_TC_ALBEDO_P223_O056884_F05_0011.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'AlbedoLocal'

    hdf = SD(FILE_NAME, SDC.READ)

    # Read dataset.
    data4D = hdf.select(DATAFIELD_NAME)

    # Convert 4-D data to 2-D data by subsetting.
    SOMBlockDim = 50
    NBandDim = 0
    data = data4D[SOMBlockDim, :, :, NBandDim].astype(np.double)

    # Read geolocation dataset from HDF-EOS2 dumper output.
    GEO_FILE_NAME = 'lat_MISR_TC_ALBEDO_P223_F05_lvl50.output'
    GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                 GEO_FILE_NAME)
    lat = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
    lat = lat.reshape(data.shape)

    GEO_FILE_NAME = 'lon_MISR_TC_ALBEDO_P223_F05_lvl50.output'
    GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                 GEO_FILE_NAME)
    lon = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
    lon = lon.reshape(data.shape)

    # Read attributes.
    attrs = data4D.attributes(full=1)
    _FillValue = attrs["_FillValue"][0]

    # Apply the fill value.
    data[data == _FillValue] = np.nan
    datam = np.ma.masked_array(data, mask=np.isnan(data))

    # Set the limit for the plot.
    m = Basemap(projection='cyl', resolution='h',
                llcrnrlat=np.min(lat), urcrnrlat=np.max(lat),
                llcrnrlon=np.min(lon), urcrnrlon=np.max(lon))
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(np.floor(np.min(lat)), np.ceil(np.max(lat)), 1),
                    labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(np.floor(np.min(lon)), np.ceil(np.max(lon)), 1),
                    labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label('No Unit')

    basename = os.path.basename(FILE_NAME)
    title = '{0}\n{1} at SOMBlockDim=50 NBandDim=0'
    plt.title(title.format(basename, DATAFIELD_NAME))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.l50.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
