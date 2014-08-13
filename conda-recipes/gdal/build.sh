#!/bin/bash

bash configure --with-python --with-hdf=$PREFIX --with-hfd5=$PREFIX --prefix=$PREFIX
make
make install

rm -rf $PREFIX/share
