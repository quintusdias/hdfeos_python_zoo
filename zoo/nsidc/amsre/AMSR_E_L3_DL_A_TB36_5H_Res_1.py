"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an NSIDC AMSR grid
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python AMSR_E_L3_DL_A_TB36_5H_Res_1.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_PYHDFEOS = True

def run(FILE_NAME):
    
    # Identify the data field.
    DATAFIELD_NAME = 'A_TB36.5H (Res 1)'

    if USE_PYHDFEOS:

        from pyhdfeos import GridFile

        gdf = GridFile(FILE_NAME)

        grid = gdf.grids['Ascending_Land_Grid']
        lat, lon = grid[:]

        field = grid.fields[DATAFIELD_NAME]
        data = field[:].astype(np.float64)
        _FillValue = field.attrs['_FillValue']

    else:
        from pyhdf.SD import SD, SDC
        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:,:].astype(np.float64)

	    # Read geolocation dataset from HDF-EOS2 dumper output.
        GEO_FILE_NAME = 'lat_AMSR_E_L3_DailyLand_V06_20050118_Ascending_Land_Grid.output'
        try:
            GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], 
                                         GEO_FILE_NAME)
        except KeyError:
            pass

        lat = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
        GEO_FILE_NAME = 'lon_AMSR_E_L3_DailyLand_V06_20050118_Ascending_Land_Grid.output'
        try: 
            GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], 
                                         GEO_FILE_NAME)
        except KeyError:
            pass
        lon = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])

        # Read attributes.
        attrs = data2D.attributes(full=1)
        fva=attrs["_FillValue"]
        _FillValue = fva[0]        


    # Apply the attributes information.
    # Ref:  http://nsidc.org/data/docs/daac/ae_land3_l3_soil_moisture/data.html
    data[data == _FillValue] = np.nan
    data *= 0.1
    data = np.ma.masked_array(data, np.isnan(data))
    long_name = DATAFIELD_NAME
    units = 'Kelvin'

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, llcrnrlon=-180, urcrnrlat=90, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, data, latlon=True)
    cb = m.colorbar()
    cb.set_label(units)
    
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)



if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'AMSR_E_L3_DailyLand_V06_20050118.hdf'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
