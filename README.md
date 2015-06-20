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

Anaconda
========
Anaconda is ideal for the Windows platform although the netcdf4
library was not compiled with hdf4 support. For HDF4, you can use
PyHDF at http://hdfeos.org/software/pyhdf.php. The H5PY package
that is installed by default will read HDF5 swath files, and the
gdal package (not installed by default) will read both HDF4 and
HDF5 grid files.

    $ conda install basemap
    $ conda install netcdf4
    $ conda install gdal

Mac
===

If you use MacPorts, you should install the hdf4 and dap variant of the netcdf
port.

    $ sudo port install netcdf +hdf4 +dap

Python2
-------
Also with MacPorts, you should install the hdf4, hdf5, and netcdf
variant of the gdal port along with basemap.

    $ sudo port install gdal +hdf4 +hdf5 +netcdf
    $ sudo port install py27-matplotlib-basemap
    $ sudo port install py27-netcdf4

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
