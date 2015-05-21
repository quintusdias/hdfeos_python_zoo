"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize an LaRC MISR SOM
grid file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MISR_ELLIPSOID_P117_F03_BlueRR_lvl0_179.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

References:

[1] http://eosweb.larc.nasa.gov/PRODOCS/misr/DPS/DPS_v50_RevS.pdf
[2] http://eospso.gsfc.nasa.gov/eos_homepae/for_scientists/atbd/docs/MISR/atbd-misr-01.pdf
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np
from pyhdf.HDF import *
from pyhdf.SD import *
from pyhdf.V import *


def run():

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    FILE_NAME = 'MISR_AM1_GRP_ELLIPSOID_GM_P117_O058421_BA_F03_0024.hdf'
    if 'HDFEOS_ZOO_DIR' in os.environ.keys():
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], FILE_NAME)

    # Identify the data field.
    DATAFIELD_NAME = 'Blue Radiance/RDQI'

    hdf = SD(FILE_NAME, SDC.READ)

    # Read dataset.
    data3D = hdf.select(DATAFIELD_NAME)
    data = data3D[:, :, :]

    # Read geolocation dataset from HDF-EOS2 dumper output.
    GEO_FILE_NAME = 'lat_MISR_ELLIPSOID_P117_F03.output'
    GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                 GEO_FILE_NAME)
    lat = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
    lat = lat.reshape(data.shape)

    GEO_FILE_NAME = 'lon_MISR_ELLIPSOID_P117_F03.output'
    GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'],
                                 GEO_FILE_NAME)
    lon = np.genfromtxt(GEO_FILE_NAME, delimiter=',', usecols=[0])
    lon = lon.reshape(data.shape)

    # Read attributes.
    attrs = data3D.attributes(full=1)
    _FillValue = attrs["_FillValue"][0]

    # Read scale factor attribute.
    f = HDF(FILE_NAME, HC.READ)
    v = f.vgstart()
    vg = v.attach(8)

    # PyHDF cannot read attributes from Vgroup properly.
    # sfa = vg.attr('Scale Factor')
    # scale_factor = sfa.get()
    vg.detach()
    v.end()

    # Set it manually using HDFView.
    scale_factor = 0.047203224152326584

    # We need to shift bits for "RDQI" to get "Blue Band "only.
    # See the page 84 of "MISR Data Products Specifications (rev. S)".
    # The document is available at [1].
    datas = np.right_shift(data, 2)
    dataf = datas.astype(np.double)

    # Apply the fill value.
    dataf[data == _FillValue] = np.nan

    # Filter out values (> 16376) used for "Flag Data".
    # See Table 1.2 in "Level 1 Radiance Scaling and Conditioning
    # Algorithm  Theoretical Basis" document [2].
    dataf[datas > 16376] = np.nan
    datam = np.ma.masked_array(dataf, mask=np.isnan(dataf))

    # Apply scale facotr.
    datam = scale_factor * datam

    nblocks = data.shape[0]
    ydimsize = data.shape[1]
    xdimsize = data.shape[2]

    datam = datam.reshape(nblocks*ydimsize, xdimsize)
    lat = lat.reshape(nblocks*ydimsize, xdimsize)
    lon = lon.reshape(nblocks*ydimsize, xdimsize)

    # Set the limit for the plot.
    m = Basemap(projection='cyl', resolution='h',
                llcrnrlat=np.min(lat), urcrnrlat=np.max(lat),
                llcrnrlon=np.min(lon), urcrnrlon=np.max(lon))
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 120, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 45), labels=[0, 0, 0, 1])
    m.pcolormesh(lon, lat, datam, latlon=True)
    cb = m.colorbar()
    cb.set_label(r'$Wm^{-2}sr^{-1}{\mu}m^{-1}$')

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, 'Blue Radiance'))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":
    run()
