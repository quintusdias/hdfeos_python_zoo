"""
This example code illustrates how to access and visualize a GESDISC MLS file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MLS_L2GP_v01_L2gpValue.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

References
[1] http://mls.jpl.nasa.gov/data/v2-2_data_quality_document.pdf
"""

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
        data_grp = dset.groups['HDFEOS'].groups['SWATHS'].groups['BrO'].groups['Data Fields']
        data_var = data_grp.variables['L2gpValue']

        data = data_var[399,:]
        units = data_var.Units
        fill_value = data_var._FillValue
        missing_value = data_var.MissingValue
        title = data_var.Title

        data[data == fill_value] = np.nan
        data[data == missing_value] = np.nan
        data = np.ma.masked_array(data, np.isnan(data))

        geo_grp = dset.groups['HDFEOS'].groups['SWATHS'].groups['BrO'].groups['Geolocation Fields']
        pres_var = geo_grp.variables['Pressure']
        pressure = pres_var[:]
        pres_units = pres_var.Units

        time_var = geo_grp.variables['Time']
        time = time_var[:]

    else:
        
        import h5py

        path = '/HDFEOS/SWATHS/BrO'
        DATAFIELD_NAME = path + 'CloudFraction'
        with h5py.File(FILE_NAME, mode='r') as f:
            varname = path + '/Data Fields/BrO'
            dset = f[varname]
            data = dset[399, :]

            # Retrieve any attributes that may be needed later.
            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            missing_value = f[varname].attrs['MissingValue']
            fill_value = f[varname].attrs['_FillValue']
            title = f[varname].attrs['Title'].decode()
            units = f[varname].attrs['Units'].decode()

            # Retrieve the geolocation data.
            varname = path + '/Geolocation Fields/Pressure'
            pressure = f[varname][:]
            pres_units = f[varname].attrs['Units'].decode()

            varname = path + '/Geolocation Fields/Time'
            time = f[varname][:]

    # Convert to minutes, time from start.
    time = (time - time[0]) / 60.0

    # Create an "elapsed time" variable (International Atomic Time).
    time1lvl = str(int(time[399]))

    # Read MLS Data Quality Document [1] for useful range in BrO data, which is
    # 3.2hPa - 10hPa
    plt.plot(pressure[12:16], data[12:16])

    plt.xlabel('Pressure ({0})'.format(pres_units))
    plt.ylabel('{0} ({1})'.format(title, units))
    fig = plt.gcf()
    
    plt.title('{0} at Time = {1}'.format(title, time1lvl))
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, 'L2gpValue')
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MLS-Aura_L2GP-BrO_v02-23-c01_2010d255.he5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
