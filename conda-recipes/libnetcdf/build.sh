#!/bin/bash

export CFLAGS="-I$PREFIX/include $CFLAGS"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS -lmfhdf -ldf -ljpeg"

./configure \
    --enable-shared \
    --enable-hdf4 \
    --enable-hdf4-file-tests \
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
