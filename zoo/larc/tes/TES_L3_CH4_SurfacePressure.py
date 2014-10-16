"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a TES L3 Grid file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python TES_L3_CH4_SurfacePressure.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os
import re

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

def run(FILE_NAME):
    
    with h5py.File(FILE_NAME, mode='r') as f:

        # Need to retrieve the grid metadata.  The hdfeos5 library stores it
        # in a string dataset.
        METADATA_FIELD = '/HDFEOS INFORMATION/StructMetadata.0'
        gridmeta = str(f[METADATA_FIELD][...])

        # Need to transpose the data
        DATA_FIELD = '/HDFEOS/GRIDS/NadirGrid/Data Fields/SurfacePressure'
        data = f[DATA_FIELD][...].astype(np.float64).T
        fillvalue = f[DATA_FIELD].attrs['_FillValue']
        missingvalue = f[DATA_FIELD].attrs['MissingValue']
        title = f[DATA_FIELD].attrs['Title']
        units = f[DATA_FIELD].attrs['Units']

    invalid = np.logical_or(data == fillvalue[0], data == missingvalue[0])
    data[invalid] = np.nan
    data = np.ma.masked_array(data, np.isnan(data))

    # Construct the grid.  Use regular expressions to tease out the
    # extents of the grid.  
    ul_regex = re.compile(r'''UpperLeftPointMtrs=\(
                              (?P<upper_left_x>[+-]?\d+\.\d+)
                              ,
                              (?P<upper_left_y>[+-]?\d+\.\d+)
                              \)''', re.VERBOSE)
    match = ul_regex.search(gridmeta)
    x0 = np.float(match.group('upper_left_x')) / 1e6
    y0 = np.float(match.group('upper_left_y')) / 1e6

    lr_regex = re.compile(r'''LowerRightMtrs=\(
                              (?P<lower_right_x>[+-]?\d+\.\d+)
                              ,
                              (?P<lower_right_y>[+-]?\d+\.\d+)
                              \)''', re.VERBOSE)
    match = lr_regex.search(gridmeta)
    x1 = np.float(match.group('lower_right_x')) / 1e6
    y1 = np.float(match.group('lower_right_y')) / 1e6
        
    ny, nx = data.shape
    x = np.linspace(x0, x1, nx)
    y = np.linspace(y0, y1, ny)
    lon, lat = np.meshgrid(x, y)
    
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-45, 91, 45), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180, 45), labels=[0, 0, 0, 1])

    m.pcolormesh(lon, lat, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)


    fig = plt.gcf()
    # plt.show()
    
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, title))

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'TES-Aura_L3-CH4_r0000010410_F01_07.he5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    

