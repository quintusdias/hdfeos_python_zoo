"""
This example code illustrates how to access and visualize a GESDISC GOSAT ACOS
L2 Swath HDF5 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python acos_L2s_110101_02_Production_v110110_L2s2800_r01_PolB.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

You can get the meaningful description of data from the "README Document for
ACOS Level 2 Standard Product" [1] or the XML description file that isprovided
by GES-DISC.

References
----------
[1] ftp://aurapar1u.ecs.nasa.gov/ftp/data/s4pa/GOSAT_TANSO_Level2/ACOS_L2S.002/doc/README.ACOS_L2S_v2.8.pdf
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
        data_var = dset.groups['RetrievalResults'].variables['xco2']
        lat_var = dset.groups['SoundingGeometry'].variables['sounding_latitude_geoid']
        lon_var = dset.groups['SoundingGeometry'].variables['sounding_longitude_geoid']

        # Read the data.
        #data = data_var[0,:,:]
        lat = lat_var[:]
        lon = lon_var[:]
        #lev = lev_var[:]
        #date = date_var[0]

        # Read the needed attributes.
        #data_units = data_var.units
        #lat_units = lat_var.units
        #lev_units = lev_var.units
        #data_longname = data_var.long_name
        #lat_longname = lat_var.long_name
        #lev_longname = lev_var.long_name

    else:
        
        import h5py

        with h5py.File(FILE_NAME, mode='r') as f:
            dset_var = f['/RetrievalResults/xco2']
            dset_lat = f['/SoundingGeometry/sounding_latitude_geoid']
            dset_lon = f['/SoundingGeometry/sounding_longitude_geoid']
            dset_lev = f['/SoundingGeometry/sounding_altitude']
            dset_time = f['/SoundingHeader/sounding_time_tai93']

            # Read the data.
            #data = dset_var[0,:,:]
            lat = dset_lat[:]
            lon = dset_lon[:]
            #lev = dset_lev[:]
            #date = dset_date[0]

            # Read the needed attributes.
            #data_units = dset_var.attrs['units']
            #lat_units = dset_lat.attrs['units']
            #lev_units = dset_lev.attrs['units']
            #data_longname = dset_var.attrs['long_name']
            #lat_longname = dset_lat.attrs['long_name']
            #lev_longname = dset_lev.attrs['long_name']

            # H5PY doesn't automatically turn the data into a masked array.
            #fillvalue = dset_var.attrs['_FillValue']
            #data[data == fillvalue] = np.nan
            #data = np.ma.masked_array(data, np.isnan(data))

    # Draw an orthographic projection using the low resolution coastline
    # database.
    m = Basemap(projection='ortho', resolution='l', lat_0=-55, lon_0 = 120)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-80., -0., 20.))
    m.drawmeridians(np.arange(-180., 181., 20.))
    x, y = m(lon, lat)
    m.plot(x, y)

    # Annotate the starting point.
    m.plot(x[0], y[0], marker='o', color='red')
    plt.annotate('START',
                 xy=(x[0] + 195000, y[0] + 5000),
                 xycoords='data',
                 color='red')

    fig = plt.gcf()
    
    plt.title('Trajectory')
    plt.show()
    
    png = "{0}.{1}.png".format(os.path.basename(FILE_NAME)[:-4], 'trajectory')
    fig.savefig(png)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'acos_L2s_110101_02_Production_v110110_L2s2800_r01_PolB_110124184213.h5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
