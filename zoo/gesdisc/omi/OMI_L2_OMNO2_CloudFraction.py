"""
This example code illustrates how to access and visualize a GESDISC OMI file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python OMI_L2_OMNO2_CloudFraction.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = True

def run(FILE_NAME):
    DATAFIELD_NAME = 'CloudFraction'
    
    if USE_NETCDF4:

        from netCDF4 import Dataset

        dset = Dataset(FILE_NAME)
        grp = dset.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2']
        var = grp.groups['Data Fields'].variables[DATAFIELD_NAME]
        
        # netCDF4 doesn't quite handle the scaling correctly in this case since
        # the attributes are "ScaleFactor" and "AddOffset" instead of
        # "scale_factor" and "add_offset".  We'll do the scaling and conversion
        # to a masked array ourselves.
        var.set_auto_maskandscale(False)

        data = var[:].astype(np.float64)
    
        # Retrieve any attributes that may be needed later.
        scale = var.ScaleFactor
        offset = var.Offset
        title = var.Title
        missing_value = var.MissingValue
        fill_value = var._FillValue

        # Retrieve the geolocation data.
        latitude = grp.groups['Geolocation Fields'].variables['Latitude'][:]
        longitude = grp.groups['Geolocation Fields'].variables['Longitude'][:]

    else:
        
        import h5py

        path = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/'
        DATAFIELD_NAME = path + 'CloudFraction'
        with h5py.File(FILE_NAME, mode='r') as f:
            dset = f[DATAFIELD_NAME]
            data =dset[:].astype(np.float64)

            # Retrieve any attributes that may be needed later.
            scale = f[DATAFIELD_NAME].attrs['ScaleFactor']
            offset = f[DATAFIELD_NAME].attrs['Offset']
            missing_value = f[DATAFIELD_NAME].attrs['MissingValue']
            fill_value = f[DATAFIELD_NAME].attrs['_FillValue']
            title = f[DATAFIELD_NAME].attrs['Title']

            # Retrieve the geolocation data.
            path = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/'
            latitude = f[path + 'Latitude'][:]
            longitude = f[path + 'Longitude'][:]

    data[data == missing_value] = np.nan
    data[data == fill_value] = np.nan
    data = scale * (data - offset)
    datam = np.ma.masked_where(np.isnan(data), data)

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    
    # Render the image in the projected coordinate system.
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, datam)
    m.colorbar()
    fig = plt.gcf()
    
    plt.title('{0})'.format(title))
    plt.show()
    
    png = "{0}.{1}.png".format(os.path.basename(FILE_NAME)[:-4],
                               os.path.basename(DATAFIELD_NAME))
    fig.savefig(png)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'OMI-Aura_L2-OMNO2_2008m0720t2016-o21357_v003-2008m0721t101450.he5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
