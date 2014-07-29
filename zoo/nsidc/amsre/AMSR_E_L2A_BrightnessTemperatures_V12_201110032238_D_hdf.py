"""
This example code illustrates how to access and visualize an NSIDC AMSR_E 
version 3 L2A HDF-EOS2 swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AMSR_E_L2A_BrightnessTemperatures_V12_201110032238_D_hdf.py

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

    DATAFIELD_NAME = '89.0V_Res.5B_TB_(not-resampled)'
    
    nc = Dataset(FILE_NAME)
    data = nc.variables[DATAFIELD_NAME][:].astype(np.float64)
    
    # Replace the filled value with NaN, replace with a masked array.
    # Apply the scaling equation.  These attributes are named in a VERY
    # non-standard manner.
    scale_factor = getattr(nc.variables[DATAFIELD_NAME], 'SCALE FACTOR')
    add_offset = nc.variables[DATAFIELD_NAME].OFFSET

    data[data == -32768] = np.nan
    data = data * scale_factor + add_offset
    datam = np.ma.masked_array(data, np.isnan(data))
    
    # There is a wrap-around effect to deal with, as some of the swath extends
    # eastward over the international dateline.  Adjust the longitude to avoid
    # the swath being smeared.
    latitude = nc.variables['Latitude'][:]
    longitude = nc.variables['Longitude'][:]
    longitude[longitude < -150] += 360

    # Render the plot in a global projection.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-150, urcrnrlon = 210)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 30), [1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181., 45), [0, 0, 0, 1])
    x, y = m(longitude, latitude)
    m.pcolormesh(x, y, datam)
    m.colorbar()
    
    units = 'Kelvin'
    plt.title('{0} ({1})'.format(DATAFIELD_NAME, units))

    fig = plt.gcf()
    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'AMSR_E_L2A_BrightnessTemperatures_V12_201110032238_D.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
    
