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

    python CER_ZAVG_Aqua_1_0_Z_MHA_CLH_IPD_lvl0.py

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
    data = data4D[0,:,:,0].astype(np.float64)

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
    
    # Generate Monthly Hourly Avgs which is 3 hour interval (3*8hr = 24hrs) [2].
    MHA = np.linspace(1, 8, 8)

    plt.contourf(MHA, lat, datam.T)
    plt.ylabel('Latitude (degrees_north)')
    plt.xlabel('Monthly 3-hourly')
    plt.xticks(MHA, ['00-03 GMT','03-06 GMT','06-09 GMT','09-12 GMT','12-15 GMT','15-18 GMT','18-21 GMT','21-24 GMT'], fontsize='8')
    cb = plt.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)

    # I guess Stats=0 means "mean" and Stats=1 means "standard
    # deviation." See [3] and page 154 of [4].
    plt.title('{0}\n{1}'.format(basename, '/1.0 Degree Zonal/Monthly Hourly Averages/Cloud Layer High/\n Ice Particle Diameter (Mean)'), fontsize=11)
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
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
# [2] http://eosweb.larc.nasa.gov/GUIDE/dataset_documents/cer_syn-avg-zavg.html
# [3] http://eosweb.larc.nasa.gov/PRODOCS/ceres/readme/readme_cer_zavg_SampleRead_R5-691.txt
# [4] http://ceres.larc.nasa.gov/documents/DPC/DPC_current/pdfs/DPC_all.pdf
