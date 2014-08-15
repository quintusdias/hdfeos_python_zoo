#!/bin/bash

bash configure --with-python \
    --with-hdf4=$PREFIX \
    --with-hdf5=$PREFIX \
    --with-geos=YES \
    --prefix=$PREFIX
make
make install

#rm -rf $PREFIX/share
mkdir -p $PREFIX/share/gdal
cp data/csv $PREFIX/share/gdal
cp data/*.wkt $PREFIX/share/gdal
