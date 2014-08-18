Anaconda
========

Mac, Linux
----------

    $ conda install conda-build binstar
    $ conda create -n hdfeos_build python numpy
    $ conda build hdf4 hdf5 h5py geos libnetcdf netcdf4-python gdal
    $ binstar upload /home/jevans/anaconda3/conda-bld/linux-64/hdf4-4.2.10-1.tar.bz2
    $ binstar upload /home2/jevans/anaconda3/conda-bld/linux-64/hdf5-1.8.13-1.tar.bz2
    $ binstar upload /home2/jevans/anaconda3/conda-bld/linux-64/h5py-2.3.1-np18py34_0.tar.bz2
    $ binstar upload /home/jevans/anaconda3/conda-bld/linux-64/libnetcdf-4.2.1.1-2.tar.bz2
    $ binstar upload /home2/jevans/anaconda3/conda-bld/linux-64/netcdf4-python-1.0.9-py34_1.tar.bz2
    $ binstar upload /home/jevans/anaconda3/conda-bld/linux-64/gdal-1.10.1-py34_1.tar.bz2
