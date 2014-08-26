"""
This example code illustrates how to access and visualize an NSIDC ICESat/GLAS
L2 swath data file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python GLAH13_633_2103_001_1317_0_01_0001.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

# Can do this using either netCDF4 or h5py.
USE_NETCDF4 = True

def run(FILE_NAME):

    DATAFIELD_NAME = '/Data_1HZ/Geolocation/d_'

    if USE_NETCDF4:
    
        from netCDF4 import Dataset
    
        nc = Dataset(FILE_NAME)

        latvar = nc.groups['Data_1HZ'].groups['Geolocation'].variables['d_lat']
        latitude = latvar[:]
        lat_vr = [latvar.valid_min, latvar.valid_max]

        lonvar = nc.groups['Data_1HZ'].groups['Geolocation'].variables['d_lon']
        longitude = lonvar[:]
        lon_vr = [lonvar.valid_min, lonvar.valid_max]

        tempvar = nc.groups['Data_1HZ'].groups['Atmosphere'].variables['d_Surface_temp']
        temp = tempvar[:]
        temp_vr = [tempvar.valid_min, tempvar.valid_max]

        timevar = nc.groups['Data_1HZ'].groups['Time'].variables['d_UTCTime_1']
        time = timevar[:]

    else:
    
        import h5py
    
        with h5py.File(FILE_NAME, mode='r') as f:
    
            latitude = f['/Data_1HZ/Geolocation/d_lat'][:]
            longitude = f['/Data_1HZ/Geolocation/d_lon'][:]
            temp = f['/Data_1HZ/Atmosphere/d_Surface_temp'][:]
            time = f['/Data_1HZ/Time/d_UTCTime_1'][:]
    
            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            long_name = f[DATAFIELD_NAME].attrs['long_name'].decode()

    # apply valid max, min
    latitude[latitude < lat_vr[0]] = np.nan
    latitude[latitude > lat_vr[1]] = np.nan
    longitude[longitude < lon_vr[0]] = np.nan
    longitude[longitude > lon_vr[1]] = np.nan
    temp[temp < temp_vr[0]] = np.nan
    temp[temp > temp_vr[1]] = np.nan

    fig = plt.figure(figsize=(15, 6))
    plt.subplot(1, 2, 1)

    # Plot surface temperature vs time.
    plt.plot(time, temp)

    # Plot the trajectory
    plt.subplot(1, 2, 2)

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(-180, 180., 45.))
    m.pcolormesh(longitude, latitude, data, latlon=True, vmin=0, vmax=1)
    m.colorbar()
    plt.title('{0}'.format(long_name))

    fig = plt.gcf()
    plt.show()

    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'GLAH13_633_2103_001_1317_0_01_0001.h5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
