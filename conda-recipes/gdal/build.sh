#!/bin/bash

bash configure --with-python --with-hdf4=$PREFIX --with-hdf5=$PREFIX --prefix=$PREFIX
make
make install

rm -rf $PREFIX/share
