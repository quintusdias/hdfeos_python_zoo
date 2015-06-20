"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LaRC CERES ES4 TRMM
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python CER_ES4_TRMM_Longwave_Flux_2_5_R.py

The netCDF file must either be in your current working directory
or in a directory specified by the environment variable HDFEOS_ZOO_DIR.

In order for the netCDF code path to work, the netcdf library must be compiled
with HDF4 support.  Please see the README for details.

References
----------
[1] http://ceres.larc.nasa.gov/documents/collect_guide/pdf/ES4_CG_R1V1.pdf
[2] http://eosweb.larc.nasa.gov/PRODOCS/ceres/images/TRMM-PFM.html
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
    FILE_NAME = 'CER_ES4_TRMM-PFM_Edition1_009001.199808.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    DATAFIELD_NAME = 'Longwave Flux (2.5R)'

    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        data = nc.variables[DATAFIELD_NAME][:].astype(np.float64)

    else:

        from pyhdf.SD import SD, SDC

        hdf = SD(FILE_NAME, SDC.READ)

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:]

    # Set fillvalue and units.
    # See "CERES Data Management System ES-4 Collection Guide" [1] and a sample
    # image by NASA [2] for details.  The fillvalue is 3.4028235E38.  Here, we
    # just use the max of the data.
    fillvalue = np.max(data)
    data[data == fillvalue] = np.nan
    datam = np.ma.masked_array(data, mask=np.isnan(data))

    # Set fillvalue and units.
    # See "CERES Data Management System ES-4 Collection Guide" [1] and a
    # sample image by NASA [2] for details.
    # The fillvalue is 3.4028235E38. Here, we use max value from the dataset.
    units = 'Watts/Meter^2'
    ysize, xsize = data.shape
    xinc = 360.0 / xsize
    yinc = 180.0 / ysize
    x0, x1 = (-180, 180)
    y0, y1 = (-90, 90)
    longitude = np.linspace(x0 + xinc/2, x1 - xinc/2, xsize)
    latitude = np.linspace(y0 + yinc/2, y1 - yinc/2, ysize)

    # Flip the latitude to run from 90 to -90
    latitude = latitude[::-1]

    # The data is global, so render in a global projection.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 90, 45))
    m.drawmeridians(np.arange(-180., 180, 45),
                    labels=[True, False, False, True])
    m.pcolormesh(longitude, latitude, datam, latlon=True)
    cb = m.colorbar()

    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, DATAFIELD_NAME))
    fig = plt.gcf()
    plt.show(block=False)
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)

if __name__ == "__main__":
    run()

# References
# [1] http://ceres.larc.nasa.gov/documents/collect_guide/pdf/ES4_CG_R1V1.pdf
