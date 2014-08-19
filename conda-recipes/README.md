Anaconda
========

Mac, Linux
----------

    $ conda install conda-build binstar
    $ conda create -n hdfeos_build python
    $ conda build hdf4 hdf5 h5py libnetcdf netcdf4-python proj4 geos gdal
    $ binstar upload /home/jevans/anaconda3/conda-bld/linux-64/hdf4-4.2.10-1.tar.bz2
    $ binstar upload /home2/jevans/anaconda3/conda-bld/linux-64/hdf5-1.8.13-1.tar.bz2
    $ binstar upload /home2/jevans/anaconda3/conda-bld/linux-64/h5py-2.3.1-np18py34_0.tar.bz2
    $ binstar upload /home/jevans/anaconda3/conda-bld/linux-64/libnetcdf-4.2.1.1-2.tar.bz2
    $ binstar upload /home2/jevans/anaconda3/conda-bld/linux-64/netcdf4-python-1.0.9-py34_1.tar.bz2
    $ binstar upload /home2/jevans/anaconda3/conda-bld/linux-64/proj4-4.8.0-0.tar.bz2
    $ binstar upload /home2/jevans/anaconda3/conda-bld/linux-64/gdal-1.11.0-np18py34_0.tar.bz2


