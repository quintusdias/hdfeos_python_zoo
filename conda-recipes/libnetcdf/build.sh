#!/bin/bash

export CFLAGS="-I$PREFIX/include $CFLAGS"
export CPPFLAGS="-I$PREFIX/include $CPPFLAGS"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS -lmfhdf -ldf -ljpeg"

ARCH="$(uname 2>/dev/null)"
case ${ARCH} in
    'Darwin')
        export DYLD_FALLBACK_LIBRARY_PATH=$PREFIX/lib
        ;;
    'Linux')
        ;;
    *)
        echo -e "Unsupported machine type: ${ARCH}";
        exit 1;
        ;;
esac


./configure \
    --enable-shared \
    --enable-netcdf-4 \
    --enable-hdf4 \
    --enable-dap \
    --without-ssl \
    --without-libidn \
    --disable-ldap \
    --prefix=$PREFIX
make
make install

rm -rf $PREFIX/share
