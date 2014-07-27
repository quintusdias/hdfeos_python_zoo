"""
This example code illustrates how to access and visualize a GESDISC TRMM file
in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python TRMM_2A12_cldWater_lvl9.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):

    DATAFIELD_NAME = 'cldWater'
    
    nc = Dataset(FILE_NAME)
    var = nc.variables[DATAFIELD_NAME]

    # cldWater has "scale_factor" and "add_offset" attributes, but the scaling
    # equation to be used here is not 
    #
    #     data = data * scale + offset
    #
    # We'll turn autoscaling off in order to correctly scale the data.  Also,
    # the fill value is not explicitly set, but appears to be -9999.
    var.set_auto_maskandscale(False)
    data = var[:,:,9].astype(np.double)
    data[data == -9999] = np.nan
    data = data / var.scale_factor + var.add_offset
    datam = np.ma.masked_array(data, np.isnan(data))
    
    # Retrieve the geolocation data.
    latitude = nc.variables['geolocation'][:,:,0]
    longitude = nc.variables['geolocation'][:,:,1]
    
    # There is a wrap-around effect to deal with.  Adjust the longitude by
    # modulus 360 to avoid the swath being smeared.
    longitude[longitude < 0] += 360
    longitude[longitude > 310] -= 360
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-50, urcrnrlon = 310)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-45, 315., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, datam, latlon=True, vmin=0, vmax=0.192)
    m.colorbar()
    plt.title('{0} (g/m^3)'.format(DATAFIELD_NAME))

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = '2A12.100402.70512.6.HDF'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
