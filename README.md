============
Instructions
============
The primary method used in the examples for reading HDF-EOS2 files (HDF4 format)
is to use the netCDF4 package.  Most of the time, however, this
package does not come ready-out-of-the-box to read HDF4 files because
the underlying netcdf library does not compile with hdf4 support.  The
directions below specify how this can be overcome in certain situations.

In cases where the datafile is HDF-EOS5, HDF4 support is not needed.  The
netCDF4 package can usually read these files, but code is also provided for
reading the file with h5py.

GDAL is used to read some HDF-EOS grid files (both version 2 and 5).

Fedora 20
=========
The Fedora 20 netcdf RPM is built with hdf4 support, so therefore netcdf4-python
can read HDF4 files out of the box.  This is the ideal situation, yay Fedora!

The Fedora 20 repositories do not include a Python3 RPM for basemap, but you can
still install basemap via pip::

    $ pip-python3 install basemap --user

The RPM for the GDAL library includes support for HDF-EOS2.

Anaconda
========

Mac, Linux
----------
The following steps worked with Anaconda 2.0.1 (Python2 and 3) on Fedora 19,
CentOS 6.5, and Mac.

Download and install the free version of anaconda, then follow these
directions ::

    $ conda install basemap gdal netcdf4
    $ conda create -n hdfeos ipython matplotlib h5py basemap netcdf4 gdal
    $ conda install -c https://conda.binstar.org/jevans hdf4=4.2.10
    $ conda remove libnetcdf
    $ conda install -c https://conda.binstar.org/jevans libnetcdf=4.2.1.1

You may have to deal with issue#32 as described at
https://github.com/ContinuumIO/anaconda-issues/issues/32.  This may or
may not occur on other platforms.  For instance, it did not occur on a
CenOS 6.5 platform.

Windows
-------
Anaconda is ideal for the Windows platform although the netcdf4 library
was not compiled with hdf4 support.  The H5PY package that is installed by
default will read HDF5 swath files, and the gdal package (not installed
by default) will read both HDF4 and HDF5 grid files.  The HDF4 swath
files cannot currently be read.

Mac
===
If you use MacPorts, you should install the hdf4 and dap variant of the netcdf
port and the hdf4, hdf5, and netcdf variant of the gdal port::

    $ sudo port install netcdf +hdf4 +dap
    $ sudo port install gdal +hdf4 +hdf5 +netcdf


