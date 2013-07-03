#!/bin/bash

# Check that the user has a valid LIBMESH_DIR set
if [ -z "$LIBMESH_DIR" ]; then
  echo "LIBMESH_DIR needs to be set to the location of your libmesh checkout for this script to work!"
  exit 1;
fi

cd $LIBMESH_DIR

# Note: even though we won't run 'make install' still specify a prefix!
./configure --prefix=${LIBMESH_DIR} \
            --with-methods="${METHODS}" \
            --enable-legacy-using-namespace \
            --enable-openmp \
            --disable-tbb \
            --disable-blocked-storage
#./configure --prefix=${LIBMESH_DIR} --with-methods="${METHODS}" --enable-legacy-include-paths --enable-legacy-using-namespace  --enable-netcdf=new --enable-exodus=new --enable-nemesis=new --disable-hdf5

make -j $JOBS

# "Fake" make install - do one extra step which (so far!) lets us avoid a full-blown make install
#mkdir -p lib

#for i in $METHODS
#do
#  ${LIBMESH_DIR}/libtool --mode=install install -c libmesh_$i.la ${LIBMESH_DIR}/lib
#done

#${LIBMESH_DIR}/libtool --mode=install install -c ${LIBMESH_DIR}/contrib/netcdf/4.2.1.1/liblib/libnetcdf.la ${LIBMESH_DIR}/lib

./install_it.sh
