mkdir lib

for i in $METHODS
do
  ${LIBMESH_DIR}/libtool --mode=install install -c libmesh_$i.la ${LIBMESH_DIR}/lib
done

${LIBMESH_DIR}/libtool --mode=install install -c ${LIBMESH_DIR}/contrib/netcdf/4.3.0/liblib/libnetcdf.la ${LIBMESH_DIR}/lib

ln -s ${LIBMESH_DIR}/include/libmesh/netcdf.h ${LIBMESH_DIR}/include/netcdf.h