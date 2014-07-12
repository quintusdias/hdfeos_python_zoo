"""
This example code illustrates how to access and visualize a GESDISC HIRDLES
Zonal Average HDF-EOS5 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python HIRDLS_Aura_L3ZAD_v06_00_00_c02_2005d022_2008d077.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import datetime
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = False

def run(FILE_NAME):
    
    if USE_NETCDF4:

        from netCDF4 import Dataset

        dset = Dataset(FILE_NAME)
        grp = dset.groups['HDFEOS'].groups['ZAS'].groups['HIRDLS'].groups['Data Fields']
        data_var = grp.variables['NO2Day']
        lat_var = grp.variables['Latitude']
        lev_var = grp.variables['Pressure']
        date_var = grp.variables['Time']

        # Read the data.
        data = data_var[0,:,:]
        lat = lat_var[:]
        lev = lev_var[:]
        time = date_var[0]

        # Read the needed attributes.
        lat_units = lat_var.Units
        lev_units = lev_var.Units
        data_title = data_var.Title
        lat_title = lat_var.Title
        lev_title = lev_var.Title

    else:
        
        import h5py

        with h5py.File(FILE_NAME, mode='r') as f:
            dset_var = f['/HDFEOS/ZAS/HIRDLS/Data Fields/NO2Day']
            dset_lat = f['/HDFEOS/ZAS/HIRDLS/Data Fields/Latitude']
            dset_lev = f['/HDFEOS/ZAS/HIRDLS/Data Fields/Pressure']
            dset_date = f['/HDFEOS/ZAS/HIRDLS/Data Fields/Time']

            # Read the data.
            data = dset_var[0,:,:]
            lat = dset_lat[:]
            lev = dset_lev[:]
            time = dset_date[0]

            # Read the needed attributes.
            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            lat_units = dset_lat.attrs['Units'].decode()
            lev_units = dset_lev.attrs['Units'].decode()
            data_title = dset_var.attrs['Title'].decode()
            lat_title = dset_lat.attrs['Title'].decode()
            lev_title = dset_lev.attrs['Title'].decode()

            # H5PY doesn't automatically turn the data into a masked array.
            fillvalue = dset_var.attrs['_FillValue']
            data[data == fillvalue] = np.nan
            data = np.ma.masked_array(data, np.isnan(data))

    # The date is stored as a six-digit number, YYYYMM.  Convert it into
    # a string.
    datestr = datetime.datetime(1993,1,1) + datetime.timedelta(seconds=time)

    # Apply log scale along the y-axis to get a better image.
    lev = np.log10(lev)
    plt.contourf(lat, lev, data)
    plt.colorbar()

    plt.xlabel('{0} ({1})'.format(lat_title, lat_units))
    plt.ylabel('{0} ({1})\nin log10 scale'.format(lev_title, lev_units))
    fig = plt.gcf()
    
    plt.title('{0}\nDate:  {1}'.format(data_title,
        datestr.strftime('%Y-%m-%d %H:%M:%S')))
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, 'NO2Day')
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'HIRDLS-Aura_L3ZAD_v06-00-00-c02_2005d022-2008d077.he5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)

