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
---
If you have root access to your Mac and you are comfortable with the command
line, I would recommend using MacPorts instead of Anaconda since MacPorts
already has packages for HDF4.  If you wish to use Anaconda, however, follow
these instructions. 

    $ conda create -n hdfeos python=3.4
    $ source activate hdfeos
    $ conda install basemap
    $ conda install netcdf4
    $ conda install -c jevans hdf4
    $ git clone https://github.com/fhs/python-hdf4.git && cd python-hdf4
    $ export INCLUDE_DIRS=$HOME/anaconda3/envs/hdfeos/include
    $ export LIBRARY_DIRS=$HOME/anaconda3/envs/hdfeos/lib
    $ python setup.py install
    # Mac only
    $ export DYLD_FALLBACK_LIBRARY_PATH=$HOME/anaconda3/envs/hdfeos/lib

Linux
-----
If you only have access to a rather conservative Linux distribution such as 
RedHat/Centos or your linux box is just plain old, then Anaconda might be your
best bet.  Just follow the same instructions as for the Mac above, except you
do not need to set DYLD_FALLBACK_LIBRARY_PATH.  If you have a fairly cutting
edge Linux distribution such as Fedora, you might be better off going with the
system Python (you will probably still need to install python-hdf4 via git,
though).

MacPorts
========

If you use MacPorts, you should install the hdf4 and dap variant of the netcdf
port and the hdf4/hdf5/netcdf variant of the gdal port.

    $ sudo port install netcdf +hdf4 +dap
    $ sudo port install gdal +hdf4 +hdf5 +netcdf

The python modules for interfacing with the gdal and netcdf libraries should
then be able to read a broader range of file formats.
