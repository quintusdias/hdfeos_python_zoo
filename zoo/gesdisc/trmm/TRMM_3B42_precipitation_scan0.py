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

    python TRMM_3B42_precipitation_scan0.py

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

    DATAFIELD_NAME = 'precipitation'
    
    # Ignore the leading singleton dimension.
    dset = Dataset(FILE_NAME)
    data = dset.variables[DATAFIELD_NAME][0,:,:].astype(np.float64)
    
    # Consider 0.0 to be the fill value.
    # Must create a masked array where nan is involved.
    data[data == 0.0] = np.nan
    datam = np.ma.masked_where(np.isnan(data), data)
    
    
    # The lat and lon should be calculated manually.
    # More information can be found at:
    # http://disc.sci.gsfc.nasa.gov/additional/faq/precipitation_faq.shtml#lat_lon
    latitude = np.arange(-49.875, 49.875, 0.249375)
    longitude = np.arange(-179.875, 179.876, 0.25)
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    
    m.drawcoastlines(linewidth=0.5)
    
    m.drawparallels(np.arange(-90, 90, 30))
    m.drawmeridians(np.arange(-180, 180, 45))
    
    # Render the image in the projected coordinate system.
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, datam.T)
    m.colorbar()
    fig = plt.gcf()
    
    # Turn on tex for the units.
    tex_not_set = False
    if not mpl.rcParams['text.usetex']:
        tex_not_set = True
        mpl.rc('text', usetex=True)

    plt.title('{0} ({1})'.format(DATAFIELD_NAME, r'$\frac{mm}{hr}$'))
    plt.show()
    
    pngfile = "{0}.{1}.png".format(os.path.basename(FILE_NAME[:-4]),
                                   DATAFIELD_NAME)
    fig.savefig(pngfile)

    # Restore original tex settings if needed.
    if tex_not_set:
        mpl.rc('text', usetex=False)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = '3B42.100331.21.6A.HDF'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
