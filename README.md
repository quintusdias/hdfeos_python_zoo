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

Fedora 20
=========
The Fedora 20 netcdf RPM is built with hdf4 support, so therefore netcdf4-python
can read HDF4 files out of the box.  This is the ideal situation, yay Fedora!

The Fedora 20 repositories do not include a Python3 RPM for basemap, but you can
still install basemap via pip::

    $ pip-python3 install basemap --user

Fedora 19
=========
The Fedora 19 netcdf RPM was not built with the --with-hdf4 option, so you must
rebuild the RPM as follows::

1. Download the netcdf SRPM with ``yumdownloader --source netcdf``
2. Install the SRPM with ``rpm -i netcdf-4.2.1.1-5.fc19.src.rpm``
3. Change directories into your rpmbuild directory, i.e. ``cd ~/rpmbuild/SPECS``
4. Edit the spec file.  Add these two lines just underneath the ``%build`` line::

    export CPPFLAGS="-I/usr/include/hdf"
    export LDFLAGS="-L/usr/lib64/hdf -lmfhdf -ldf -ljpeg"

5.  Under the ``%global configure_opts`` section, add the following two lines::

    --enable-hdf4 \\\ 
    --enable-hdf4-file-tests \\\ 

6. Rebuild the RPMs with ``rpmbuild -bb netcdf.spec``
7. Install the newly-built RPM(s).

Ubuntu 13.10
============
The versions of the netcdf and hdf packages that come with Ubuntu are a bit 
dated, so we'll detail how to do a source install here.

From a base install, use apt-get to additionally install 

    * python-dev
    * ipython
    * python-matplotlib
    * python-mpltoolkits.basemap-data
    * python-mpltoolkits.basemap

Download HDF 4.2.10 from http://www.hdfgroup.org.  Configure and install with::

    $ ./configure --prefix=/usr/local --disable-netcdf --enable-fortran=no --enable-shared --disable-static
    $ make install

Download HDF5 1.8.13, configure and install with::

    $ ./configure --prefix=/usr/local --disable-static --enable-shared
    $ make install

Download netcdf 4.3.2 from Unidata, configure and install with::

    $ ./configure --prefix=/usr/local --disable-static --enable-shared --enable-hdf4 --enable-dap
    $ make install

Download netcdf4-1.1.0 from Pypi, configure and install with::

    $ python setup.py install --user


Anaconda
========
TODO:  The netcdf library was not compiled with hdf4 support.

Mac
===
If you use MacPorts, you should install the hdf4 variant of the netcdf port::

    $ sudo port install netcdf+hdf
