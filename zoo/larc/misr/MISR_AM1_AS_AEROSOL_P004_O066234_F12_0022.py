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

    python MISR_AM1_AS_AEROSOL_P004_O066234_F12_0022.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

 The file contains SOM projection. We need to use eosdump to generate 1D
 lat and lon and then convert them to 2D lat and lon accordingly.
 For example, run command as follows to get SOM projectoin lat/lon in ASCII.

 eos2dump -c1 MISR_AM1_AS_AEROSOL_P004_O066234_F12_0022.hdf RegParamsAer all > lat_MISR_AM1_AS_AEROSOL_P004_O066234_F12_0022.output
 eos2dump -c2 MISR_AM1_AS_AEROSOL_P004_O066234_F12_0022.hdf RegParamsAer all > lon_MISR_AM1_AS_AEROSOL_P004_O066234_F12_0022.output

 To properly display the data, the latitude/longitude must be remapped.

"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np
from pyhdf.HDF import *
from pyhdf.SD import *
from pyhdf.V import *


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'MISR_AM1_AS_AEROSOL_P004_O066234_F12_0022.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'RegBestEstimateSpectralOptDepth'

    hdf = SD(FILE_NAME, SDC.READ)

    # Read dataset.
    data4D = hdf.select(DATAFIELD_NAME)

    # Subset the Blue Band. 1=Blue, 2=Green, 3=Red, 4=NIR.
    data = data4D[:, :, :, 0].astype(np.double)

    # Read geolocation dataset from HDF-EOS2 dumper output.
    GEO_FILE_NAME = 'lat_MISR_AM1_AS_AEROSOL_P004_O066234_F12_0022.output'
    GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                 GEO_FILE_NAME)
    lat = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
    lat = lat.reshape(data.shape)

    GEO_FILE_NAME = 'lon_MISR_AM1_AS_AEROSOL_P004_O066234_F12_0022.output'
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

    nblocks = data.shape[0]
    ydimsize = data.shape[1]
    xdimsize = data.shape[2]

    datam = datam.reshape(nblocks * ydimsize, xdimsize)
    lat = lat.reshape(nblocks * ydimsize, xdimsize)
    lon = lon.reshape(nblocks * ydimsize, xdimsize)

    # Set the limit for the plot.
    m = Basemap(projection='cyl', resolution='h',
                llcrnrlat=np.min(lat), urcrnrlat=np.max(lat),
                llcrnrlon=np.min(lon), urcrnrlon=np.max(lon))
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 120, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label('No Unit')

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1} - Blue Band'.format(basename, DATAFIELD_NAME))
    fig = plt.gcf()
    plt.show(block=False)
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
