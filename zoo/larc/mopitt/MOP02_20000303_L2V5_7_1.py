"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a MOPITT HDF-EOS2 
swath file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOP02_20000303_L2V5_7_1.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

from pyhdf import HDF, SD, VS, V

def run(FILE_NAME):
    
    # Initialize the SD, V, and VS interfaces.
    hdf = HDF.HDF(FILE_NAME)
    v = hdf.vgstart()
    vs = hdf.vstart()
    sd = SD.SD(FILE_NAME)

    # Navigate the HDF-EOS2 structure and extract the geolocation data.  These
    # are vdatas, so the SD interface cannot retrieve it.  This is also why
    # the netcdf interface cannot retrieve it.
    xid = vs.find('Latitude')
    latid = vs.attach(xid)
    latid.setfields('Latitude')
    nrecs, _, _, _, _ = latid.inquire()
    latitude = latid.read(nRec=nrecs)
    latid.detach()

    lonid = vs.attach(vs.find('Longitude'))
    lonid.setfields('Longitude')
    nrecs, _, _, _, _ = lonid.inquire()
    longitude = lonid.read(nRec=nrecs)
    lonid.detach()

    # Extract the data.
    name = 'Retrieval Bottom Pressure'
    sds = sd.select(name)
    data = sds[:, 0].astype(np.float64)
    sds.endaccess()

    v.end()
    vs.end()
    sd.end()

    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90, 91, 45))
    m.drawmeridians(np.arange(-180, 180, 45), labels=[True,False,False,True])
    m.scatter(longitude, latitude, c=data, s=1, cmap=plt.cm.jet,
            edgecolors=None, linewidth=0)

    cb = m.colorbar()
    cb.set_label('hPa')

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename,  name))
    fig = plt.gcf()
    # plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)


if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MOP02-20000303-L2V5.7.1.val.hdf'

    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
