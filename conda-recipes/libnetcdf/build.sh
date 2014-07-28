#!/bin/bash

export CPPFLAGS="-I$PREFIX/include $CPPFLAGS"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS -lmfhdf -ldf -ljpeg"
export DYLD_LIBRARY_PATH=$PREFIX/lib

./configure \
    --enable-shared \
    --enable-hdf4 \
    --enable-netcdf-4 \
    --enable-dap \
    --without-ssl \
    --without-libidn \
    --disable-ldap \
    --disable-dap-remote-tests \
    --prefix=$PREFIX
make
make check
make install

rm -rf $PREFIX/share
