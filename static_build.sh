#!/bin/bash

export LIBMESH_ROOT=`pwd`
export LIBMESH_DIR=${LIBMESH_ROOT}/installed_static
export METHODS="opt dbg"
export JOBS=${JOBS:=1}

rm -rf build_static installed_static
mkdir build_static
cd build_static
../configure --with-methods="${METHODS}" \
             --prefix=$LIBMESH_DIR \
             --enable-legacy-include-paths \
             --enable-legacy-using-namespace \
             --disable-shared \
             --enable-static

make -j $JOBS
make install
