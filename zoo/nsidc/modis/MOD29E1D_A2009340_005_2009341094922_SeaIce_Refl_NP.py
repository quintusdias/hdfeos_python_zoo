"""
This example code illustrates how to access and visualize a NSIDC MODIS 4km
LAMAZ (Ease) Grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD29E1D_A2009340_005_2009341094922_SeaIce_Refl_NP.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""

import os
import re

import gdal
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
from netCDF4 import Dataset
import numpy as np

def run(FILE_NAME):
    
    # Identify the data field.
    GRID_NAME = 'MOD_Grid_Seaice_4km_North'
    DATAFIELD_NAME = 'Sea_Ice_by_Reflectance_NP'
    
    # Use netcdf to read the data itself.
    dset = Dataset(FILE_NAME)
    data = dset.variables[DATAFIELD_NAME][:]

    # Only use gdal to get the projection information.
    gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
                                                     GRID_NAME,
                                                     DATAFIELD_NAME)
    gdset = gdal.Open(gname)
    meta = gdset.GetMetadata()
    x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
    nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
    x = np.linspace(x0, x0 + xinc*nx, nx)
    y = np.linspace(y0, y0 + yinc*ny, ny)
    xv, yv = np.meshgrid(x, y)

    # Reproject into latlon
    lamaz = pyproj.Proj("+proj=laea +a=6371228 +lat_0=90 +lon_0=0 +units=m")
    wgs84 = pyproj.Proj("+init=EPSG:4326") 
    lon, lat= pyproj.transform(lamaz, wgs84, xv, yv)

    # Use a north polar azimuthal equal area projection.
    m = Basemap(projection='nplaea', resolution='l',
                boundinglat=20, lon_0=0)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 0, 15), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180, 45), labels=[0, 0, 0, 1])

    # Use a discretized colormap since we have only a few levels.
    # 0=missing data
    # 1=no decision
    # 11=night
    # 25=land
    # 37=inland water
    # 39=ocean
    # 50=cloud
    # 200=sea ice
    # 253=no input tile expected    
    # 254=non-production mask"
    lst = ['#727272',
           '#b7b7b7',
           '#ffff96',
           '#00ff00',
           '#232375',
           '#232375',
           '#63c6ff',
           '#ff0000',
           '#3f3f3f',
           '#000000']
    cmap = mpl.colors.ListedColormap(lst)
    bounds = [0, 1, 11, 25, 37, 39, 50, 200, 253, 254, 255]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # Render only a subset of the mesh.
    rows = slice(500, 4000, 5)
    cols = slice(500, 4000, 5)
    m.pcolormesh(lon[rows,cols], lat[rows,cols], data[rows,cols],
                 latlon=True, cmap=cmap, norm=norm)
    
    color_bar = plt.colorbar()
    color_bar.set_ticks([0.5, 5.5, 18, 31, 38, 44.5, 125, 226.5, 253.5, 254.5])
    color_bar.set_ticklabels(['missing', 'no decision', 'night', 'land',
                              'inland water', 'ocean', 'cloud', 'sea ice',
                              'no input tile\nexpected',
                              'non-production\nmask'])
    color_bar.draw_all()
    fig = plt.gcf()
    
    plt.title(DATAFIELD_NAME.replace('_',' '))
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, DATAFIELD_NAME)
    fig.savefig(pngfile)

    del gdset


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOD29E1D.A2000055.005.2006268025009.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
