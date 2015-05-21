"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a GESDISC TRMM 2A12
version 7 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python TRMM_2A12_20140308_92894_7_HDF.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = False


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = '2A12.20140308.92894.7.HDF'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'surfaceRain'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        # Retrieve the data.
        nc = Dataset(FILE_NAME)
        data = nc.variables[DATAFIELD_NAME][:].astype(np.float64)
        units = nc.variables[DATAFIELD_NAME].units

        # Retrieve the geolocation data.
        latitude = nc.variables['Latitude'][:]
        longitude = nc.variables['Longitude'][:]

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        ds = hdf.select(DATAFIELD_NAME)
        data = ds[:].astype(np.double)

        # Handle scale/osffset attributes.
        attrs = ds.attributes(full=1)
        units = attrs["units"][0]

        # Retrieve the geolocation data.
        lat = hdf.select('Latitude')
        latitude = lat[:]
        lon = hdf.select('Longitude')
        longitude = lon[:]

    # Construct an indexed version of the data.
    levels = [0.0, 0.1, 1.0, 10.0, 30.0]
    Z = np.zeros(data.shape, dtype=np.float64)
    for j in range(len(levels)-1):
        Z[np.logical_and(data >= levels[j], data < levels[j+1])] = j
    Z[data >= levels[-1]] = len(levels)

    # There is a wrap-around effect to deal with.  Adjust the longitude by
    # modulus 360 to avoid the swath being smeared.
    longitude[longitude < -165] += 360

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-165, urcrnrlon=197)

    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])

    # Use a discretized colormap since we have only five levels.
    colors = ['#0000ff', '#0088ff', '#8888ff', '#ff8888', '#ff0000']
    cmap = mpl.colors.ListedColormap(colors)
    bounds = np.arange(6)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    m.pcolormesh(longitude, latitude, Z, latlon=True, cmap=cmap, norm=norm)
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, DATAFIELD_NAME))

    # Adjust colorbar height.
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    color_bar = plt.colorbar(cax=cax)

    color_bar.set_ticks([0.5, 1.5, 2.5, 3.5, 4.5])
    color_bar.set_ticklabels(['0', '0.1', '1.0', '10', '30'])
    color_bar.set_label('Unit:'+units)

    fig = plt.gcf()
    # plt.show()

    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()
