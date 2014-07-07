"""
This example code illustrates how to access and visualize a GESDISC MEaSURES
Ozone Zonal Average HDF5 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python SBUV2_NOAA17_L2_SBUV2N17L2_2011m1231_v01_01_2012m0905t152911_h5.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import datetime
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = True

def run(FILE_NAME):
    
    if USE_NETCDF4:

        from netCDF4 import Dataset

        dset = Dataset(FILE_NAME)
        data_var = dset.groups['SCIENCE_DATA'].variables['ProfileO3Retrieved']
        lat_var = dset.groups['GEOLOCATION_DATA'].variables['Latitude']
        lev_var = dset.groups['ANCILLARY_DATA'].variables['PressureLevels']
        time_var = dset.variables['nTimes']

        # Read the data.
        data = data_var[:]
        lat = lat_var[:]
        lev = lev_var[:]
        time = time_var[:]

        # Read the needed attributes.
        data_units = data_var.units
        data_vmin = data_var.valid_min
        data_vmax = data_var.valid_max
        data_fillvalue = data_var._FillValue
        lat_units = lat_var.units
        lev_units = lev_var.units
        data_longname = data_var.long_name
        lat_longname = lat_var.long_name
        lev_longname = lev_var.long_name

    else:
        
        import h5py

        with h5py.File(FILE_NAME, mode='r') as f:

            dset_var = f['/SCIENCE_DATA/ProfileO3Retrieved']
            dset_lat = f['/GEOLOCATION_DATA/Latitude']
            dset_lev = f['/ANCILLARY_DATA/PressureLevels']
            dset_time = f['nTimes']

            # Read the data.
            data = dset_var[:,:]
            lat = dset_lat[:]
            lev = dset_lev[:]
            time = dset_time[:]

            # Read the needed attributes.
            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            data_units = dset_var.attrs['units'].decode()
            data_vmin = dset_var.attrs['valid_min']
            data_vmax = dset_var.attrs['valid_max']
            data_fillvalue = dset_var.attrs['_FillValue']
            lat_units = dset_lat.attrs['units'].decode()
            lev_units = dset_lev.attrs['units'].decode()
            data_longname = dset_var.attrs['long_name'].decode()
            lat_longname = dset_lat.attrs['long_name'].decode()
            lev_longname = dset_lev.attrs['long_name'].decode()

    # Apply the valid range.
    data[data < data_vmin] = np.nan
    data[data > data_vmax] = np.nan
    data[data == data_fillvalue] = np.nan
    data = np.ma.masked_array(data, np.isnan(data))

    # The time is stored as seconds since 1993-01-01
    # a string.
    start_time = datetime.datetime(1993,1,1) + datetime.timedelta(seconds=time[0])
    end_time = datetime.datetime(1993,1,1) + datetime.timedelta(seconds=time[70])

    # Apply log scale along the y-axis to get a better image.
    lev = np.log10(lev)
    plt.contourf(lat, lev, data.T, levels=np.arange(0,60,5))
    plt.colorbar()

    plt.xlabel('{0} ({1})'.format(lat_longname, lat_units))
    plt.ylabel('{0} ({1})\nin log10 scale'.format(lev_longname, lev_units))
    fig = plt.gcf()
    
    plt.title('{0} ({1})\n{2} - {3}'.format(data_longname,
                                           data_units,
                                           start_time,
                                           end_time))
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, 'ProfileO3Retrieved')
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'SBUV2-NOAA17_L2-SBUV2N17L2_2011m1231_v01-01-2012m0905t152911.h5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
