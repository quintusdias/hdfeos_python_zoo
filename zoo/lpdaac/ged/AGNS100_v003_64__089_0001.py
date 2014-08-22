"""
This example code illustrates how to access and visualize an LP DAAC ASTER GED
HDF5 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AGNS100_v003_64__089_0001.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

References:
    [1] https://lpdaac.usgs.gov/products/community_products_table/agns100
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

        nc = Dataset(FILE_NAME)
        var = nc.groups['Emissivity'].variables['Mean']

        # Subset for Band 10.
        data = var[1,:,:].astype(np.float64)

        # Retrieve the geolocation data.
        latitude = nc.groups['Geolocation'].variables['Latitude'][:]
        longitude = nc.groups['Geolocation'].variables['Longitude'][:]

    else:
        
        import h5py

        with h5py.File(FILE_NAME, mode='r') as f:

            dset = f['/Emissivity/Mean']
            data =dset[:].astype(np.float64)

            # Retrieve the geolocation data.
            latitude = f['/Geolocation/Latitude'][:]
            longitude = f['/Geolocation/Longitude'][:]

    # Apply the fillvalue and scaling (see [1])
    data[data == -9999] = np.nan
    data = 0.001 * data
    datam = np.ma.masked_where(np.isnan(data), data)

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='h',
                llcrnrlat=62.5,   urcrnrlat=64.5,
                llcrnrlon=-89.5,  urcrnrlon=-87.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(62, 65, 1), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-89, -87.5, 1), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    m.colorbar()
    plt.title('Mean Emissivity for Band 10')

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, 'CloudFraction')
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'AGNS100.v003.64.-089.0001.h5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)

