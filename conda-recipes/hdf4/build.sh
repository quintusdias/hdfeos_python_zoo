#!/bin/bash

mkdir -vp ${PREFIX}/bin;

ARCH="$(uname 2>/dev/null)"

#export CFLAGS="-I$PREFIX/include $CFLAGS"
export CPPFLAGS="-I${PREFIX}/include"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS"

DarwinInstallation() {

    export DYLD_FALLBACK_LIBRARY_PATH=$PREFIX/lib

    chmod +x configure;

    ./configure \
        --disable-netcdf \
        --enable-fortran=no \
        --with-jpeg=${PREFIX} \
        --disable-static \
        --without-szlib \
        --with-ssl \
        --with-zlib \
        --prefix=${PREFIX} || return 1;
    make || return 1;
    make check || return 1;
    make install || return 1;

    rm -rf ${PREFIX}/share/hdf4_examples;

    return 0;
}
LinuxInstallation() {

    export CFLAGS="-m64 -pipe -O2 -march=x86-64 -fPIC"
    chmod +x configure;

    ./configure \
        --disable-netcdf \
        --enable-fortran=no \
        --disable-static \
        --enable-linux-lfs \
        --with-ssl \
        --with-zlib \
        --prefix=${PREFIX} || return 1;
    make || return 1;
    make check || return 1;
    make install || return 1;

    rm -rf ${PREFIX}/share/hdf4_examples;

    return 0;
}

case ${ARCH} in
    'Darwin')
        DarwinInstallation || exit 1;
        ;;
    'Linux')
        LinuxInstallation || exit 1;
        ;;
    *)
        echo -e "Unsupported machine type: ${ARCH}";
        exit 1;
        ;;
esac
