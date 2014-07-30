#!/bin/bash

ARCH="$(uname 2>/dev/null)"

export CPPFLAGS="-I$PREFIX/include $CPPFLAGS"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS -lmfhdf -ldf -ljpeg"

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
