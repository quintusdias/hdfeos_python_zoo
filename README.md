============
Instructions
============
The primary method used in the examples for reading HDF-EOS2 files (HDF4 format)
is to use the netCDF4 package.  Most of the time, however, this
package does not come ready-out-of-the-box to read HDFEOS files either because
the underlying netcdf library was not compiled with hdf4 support or because
the underlying gdal library was not compiled with HDFEOS 2-or-5 support.  The
directions below specify how this can be overcome in certain situations.

In cases where the datafile is and HDF-EOS5 swath file, HDF4 support is
not needed.  The netCDF4 package can usually read these files, but code
is also provided for reading the file with h5py.

GDAL is used to read some HDF-EOS grid files (both version 2 and 5).

Fedora 20
=========
The Fedora 20 netcdf RPM is built with hdf4 support, so therefore
netcdf4-python can read HDF4 files out of the box, and the gdal RPM
includes support for HDFEOS, so for python2, everything just works.
This is the ideal situation, yay Fedora!

The Fedora 20 repositories do not include a Python3 RPM for basemap, but you can
still install basemap via pip::

    $ pip-python3 install basemap --user

The RPM for the GDAL library includes support for HDF-EOS2.

Linux Mint 17
=============
The packaged version of the netcdf library is a bit old (and is not built to 
handle HDF4 files), so we'll use the source version instead. 

* Remove libnetcdfc7 and libhdf4-0 packages if installed
* Download hdf-4.2.10 from the HDF Group website.
* Configure with ```./configure --prefix=/path/to/install --enable-netcdf=no --enable-shared --disable-static --disable-fortran``` then ```make & make install```
* Download netcdf-4.3.2 from Unidata
* configure with ```CPPFLAGS="-I/path/to/install/include" LDFLAGS="-L/path/to/install/lib  -lmfhdf -ldf -ljpeg -lz" ./configure --prefix=/path/to/install --enable-hdf4=yes --enable-netcdf-4=yes --disable-static --enable-shared``` then ```make & make install```
* Download netCDF4 from pypi (get at least version 1.1.1)
* Install with ```NETCDF4_DIR=/path/to/install python3 setup.py install --user```
* Download the gdal library version 1.11.0
* Configure with ```CPPFLAGS="-I/opt/local/include" LDFLAGS="-L/opt/local/lib" ./configure --prefix=/opt/local```
* Download the gdal python bindings from pypi, version 1.11.0
* ```mkdir -p /opt/local/lib/python2.7/site-packages```
* Install with ```python setup.py install --prefix=/opt/local```




Anaconda
========

Mac, Linux
----------
The following steps give you a netcdf library that allows you to read HDF4
files on Python3/Anaconda 2.0.1 on Fedora 19, CentOS 6.5, and Mac.

Download and install the free version of anaconda, then follow these
directions ::

    $ conda create -n hdfeos python
    $ source activate hdfeos
    $ conda install basemap 
    $ conda install --channel https://conda.binstar.org/jevans hdf4 hdf5 h5py libnetcdf netcdf4-python gdal

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
port.

    $ sudo port install netcdf +hdf4 +dap

Python2
-------
Also with MacPorts, you should install the hdf4, hdf5, and netcdf
variant of the gdal port along with basemap::

    $ sudo port install gdal +hdf4 +hdf5 +netcdf
    $ sudo port install py27-matplotlib-basemap

Python3
-------
Again with MacPorts, there is as of 2014-07-22 no basemap port for Python 3.4.
The following steps will work, though.

1.  If already installed, uninstall the py34-matplotlib port.
2.  Create a virtual environment using pyvenv-3.4 (with system site 
    packages) and activate it.
3.  Build and install matplotlib 1.3.1 from source.
4.  Download and install basemap 1.0.7.
5.  Download and install the python module for gdal from PyPi.
