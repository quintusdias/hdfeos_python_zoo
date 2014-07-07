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

    python OMI_OMCLDO2G.py

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

    DATAFIELD_NAME = 'CloudPressure'

    if USE_NETCDF4:
    
        from netCDF4 import Dataset

        dset = Dataset(FILE_NAME)
        grp = dset.groups['HDFEOS'].groups['GRIDS'].groups['CloudFractionAndPressure'].groups['Data Fields']
        var = grp.variables[DATAFIELD_NAME]
        
        # Turn off autoscaling so we can handle it uniformly for both h5py and
        # netcdf.
        var.set_auto_maskandscale(False)
        data = var[0,:,:]
        units = var.Units
        title = var.Title
        fill_value = var._FillValue
        
        # Retrieve the geolocation data.
        var = grp.variables['Longitude']
        var.set_auto_maskandscale(False)
        longitude = var[0,:,:]
        lon_fv = var._FillValue

        var = grp.variables['Latitude']
        var.set_auto_maskandscale(False)
        latitude = var[0,:,:]
        lat_fv = var._FillValue
        
    else:

        import h5py

        path = '/HDFEOS/GRIDS/CloudFractionAndPressure/Data Fields'
        with h5py.File(FILE_NAME, mode='r') as f:

            # String attributes actually come in as the bytes type and should
            # be decoded to UTF-8 (python3).
            varname = path + '/CloudPressure'
            data = f[varname][0,:,:]
            units = f[varname].attrs['Units'].decode()
            title = f[varname].attrs['Title'].decode()
            fill_value = f[varname].attrs['_FillValue'][0]
            
            # Retrieve the geolocation data.
            varname = path + '/Longitude'
            longitude = f[varname][0,:,:]
            lon_fv = f[varname].attrs['_FillValue'][0]
            varname = path + '/Latitude'
            latitude = f[varname][0,:,:]
            lat_fv = f[varname].attrs['_FillValue'][0]
            
    # The latitude and longitude grid is not complete and has a lot of fill values,
    # which is not the usual case.  Restrict the data in this case to a 1D array
    # of valid points.
    data = data[data != fill_value]
    longitude = longitude[longitude != lon_fv]
    latitude = latitude[latitude != lat_fv]
    
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
    cmap = plt.cm.jet
    sc = m.scatter(x, y, c=data, s=1, cmap=cmap, edgecolors=None,
            linewidth=0)
    m.colorbar(sc)
    fig = plt.gcf()
    
    plt.title('{0} ({1})'.format(title, units))
    plt.show()
    
    png = "{0}.{1}.png".format(os.path.basename(FILE_NAME)[:-4],
                               os.path.basename(DATAFIELD_NAME))
    fig.savefig(png)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'OMI-Aura_L2G-OMCLDO2G_2007m0129_v002-2007m0130t174603.he5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
