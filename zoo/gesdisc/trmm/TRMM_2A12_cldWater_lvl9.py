"""
This example code illustrates how to access and visualize GESDISC_TRMM file in
Python.  If you have any questions, suggestions, comments  on this example,
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

    DATAFIELD_NAME = 'cldWater'
    
    dset = Dataset(FILE_NAME)
    var = dset.variables[DATAFIELD_NAME]

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
    latitude = dset.variables['geolocation'][:,:,0]
    longitude = dset.variables['geolocation'][:,:,1]
    
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
    m.drawparallels(np.arange(-90., 90., 30.))
    m.drawmeridians(np.arange(-45, 315., 45.))
    
    # Render the image in the projected coordinate system.
    # More than 99% of the pixel values are less than 0.1.
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, datam, vmin=0, vmax=0.192)
    cb = m.colorbar()
    fig = plt.gcf()
    
    plt.title(DATAFIELD_NAME)

    # Turn on tex for the units.
    tex_not_set = False
    if not mpl.rcParams['text.usetex']:
        tex_not_set = True
        mpl.rc('text', usetex=True)

    cb.set_label(r'($\frac{g}{m^3}$)', rotation=0)
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
    hdffile = '2A12.100402.70512.6.HDF'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
