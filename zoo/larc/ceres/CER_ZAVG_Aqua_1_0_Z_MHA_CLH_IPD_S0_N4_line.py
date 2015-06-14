"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LaRC CERES ZAVG
 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python CER_ZAVG_Aqua_1_0_Z_MHA_CLH_IPD_S0_N4_line.py

The HDF4 file must either be in your current working directory
or in a directory specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from pyhdf.SD import SD, SDC
def run(FILE_NAME):

    # Identify the data field.
    DATAFIELD_NAME = 'Ice Particle Diameter'

    hdf = SD(FILE_NAME, SDC.READ)

    # Read dataset.
    # 
    # The file has many 'Ice Particle Diameter' dataset under
    # different Vgroup.
    # 
    # Use HDFView to look up ref number.
    index = hdf.reftoindex(38)
    data4D = hdf.select(index)
    data = data4D[0,3,:,0].astype(np.float64)

    # Read attributes.
    attrs = data4D.attributes(full=1)
    fva=attrs["_FillValue"]
    fillvalue = fva[0]
    ua=attrs["units"]
    units = ua[0]

    # Apply the fill value.
    data[data == fillvalue] = np.nan
    datam = np.ma.masked_array(data, mask=np.isnan(data))

    # The lat and lon should be calculated following [1].
    lat = np.linspace(89.5, -89.5, 180)

    plt.plot(lat, data)
    plt.xlabel('Latitude (degrees_north)')
    plt.ylabel('{0} ({1})'.format(DATAFIELD_NAME, units))

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, '/1.0 Degree Zonal/Monthly Hourly Averages/Cloud Layer High/\n Ice Particle Diameter (Mean; 09-12 GMT)'), fontsize=11)
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.line.png".format(basename)
    fig.savefig(pngfile)
    
    
if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'CER_ZAVG_Aqua-FM4-MODIS_Edition2B_007005.200503.hdf'
    try:
        fname = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        fname = hdffile

    run(fname)

# References
#
# [1] http://eosweb.larc.nasa.gov/PRODOCS/ceres/SRBAVG/Quality_Summaries/srbavg_ed2d/nestedgrid.html
