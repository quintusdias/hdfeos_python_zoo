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
port and the hdf4/hdf5/netcdf variant of the gdal port.

    $ sudo port install netcdf +hdf4 +dap
    $ sudo port install gdal +hdf4 +hdf5 +netcdf

The python modules for interfacing with the gdal and netcdf libraries should
then be able to read a broader range of formats.
