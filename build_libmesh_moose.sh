#!/bin/bash

export LIBMESH_ROOT=`pwd`
export LIBMESH_DIR=${LIBMESH_ROOT}/installed
export METHODS="opt dbg"
export JOBS=${JOBS:=1}

rm -rf build installed
mkdir build
cd build
../configure --with-methods="${METHODS}" \
             --prefix=$LIBMESH_DIR \
             --enable-legacy-include-paths \
             --enable-legacy-using-namespace \
             --disable-shared \
             --enable-static

make -j $JOBS
make install
