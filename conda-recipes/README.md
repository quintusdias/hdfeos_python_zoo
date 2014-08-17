Anaconda
========

Mac, Linux
----------

    $ conda install conda-build binstar
    $ conda create -n hdfeos_build python
    $ conda build hdf4 hdf5 geos libnetcdf gdal
    $ binstar upload /home/jevans/anaconda3/conda-bld/linux-64/hdf4-4.2.10-1.tar.bz2
    $ binstar upload /home/jevans/anaconda3/conda-bld/linux-64/libnetcdf-4.2.1.1-2.tar.bz2
    $ binstar upload /home/jevans/anaconda3/conda-bld/linux-64/gdal-1.10.1-py34_1.tar.bz2
