"""
This example code illustrates how to access and visualize a GESDISC OMI file
in Python.  If you have any questions, suggestions, comments  on this example,
please use the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to
see an  example of any other NASA HDF/HDF-EOS data product that is not 
listed in the HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo),
feel free to contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):
    DATAFIELD_NAME = 'CloudFraction'
    
    dset = Dataset(FILE_NAME)
    grp = dset.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2']
    var = grp.groups['Data Fields'].variables[DATAFIELD_NAME]
    data = var[:].astype(np.float64)
    
    # Scale the data appropriatedy.
    scale = var.ScaleFactor
    offset = var.Offset
    data = scale * (data - offset)
    
    # Replace the missing values with NaN.  Must create a masked array where
    # NaN is involved.
    missing_value = var.MissingValue
    data[data == missing_value] = np.nan
    fill_value = var._FillValue
    data[data == fill_value] = np.nan
    datam = np.ma.masked_where(np.isnan(data), data)
    
    # Retrieve the geolocation data.
    latitude = grp.groups['Geolocation Fields'].variables['Latitude'][:]
    longitude = grp.groups['Geolocation Fields'].variables['Longitude'][:]
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 90., 30.))
    m.drawmeridians(np.arange(-180, 180., 45.))
    
    # Render the image in the projected coordinate system.
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, datam, alpha=0.9)
    m.colorbar()
    fig = plt.gcf()
    
    plt.title('{0}\n{1})'.format(FILE_NAME, var.Title))
    plt.show()
    
    filename = "{0}.{1}.png".format(os.path.basename(FILE_NAME)[:-4],
                                    DATAFIELD_NAME)
    fig.savefig(filename)

if __name__ == "__main__":

    try:
        # If a certain environment variable is set, look there for the input
        # file.
        basename = 'OMI-Aura_L2-OMNO2_2008m0720t2016-o21357_v003-2008m0721t101450.he5'
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'], basename)
    except KeyError:
        # Look in the working directory.
        fname = basename

    run(fname)
    
