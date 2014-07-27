#!/bin/bash

export CFLAGS="-I$PREFIX/include $CFLAGS"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS"

# Only need this on mac.
export DYLD_LIBRARY_PATH=$PREFIX/lib

./configure \
    --enable-shared \
    --disable-static \
    --enable-fortran=no \
    --disable-netcdf \
    --prefix=$PREFIX
make
make check
make install

rm -rf $PREFIX/share
