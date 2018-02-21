// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// C++ includes
#include <unistd.h> // mkstemp
#include <fstream>

#include "libmesh/libmesh_config.h"

// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/dense_matrix.h"

namespace
{
using namespace libMesh;

namespace libMesh
{


//-----------------------------------------------------------------------
// LumpedMatrix members


// Constructor
template <typename T>
LumpedMatrix<T>::LumpedMatrix(const Parallel::Communicator & comm_in) :
  SparseMatrix<T>(comm_in),
  _destroy_mat_on_exit(true)
{}



// Constructor taking an existing Mat but not the responsibility
// for destroying it
template <typename T>
LumpedMatrix<T>::LumpedMatrix(Mat mat_in,
                            const Parallel::Communicator & comm_in) :
  SparseMatrix<T>(comm_in),
  _destroy_mat_on_exit(false)
{
  this->_mat = mat_in;
  this->_is_initialized = true;
}



// Destructor
template <typename T>
LumpedMatrix<T>::~LumpedMatrix()
{
  this->clear();
}


template <typename T>
void LumpedMatrix<T>::init (const numeric_index_type m_in,
                           const numeric_index_type n_in,
                           const numeric_index_type m_l,
                           const numeric_index_type n_l,
                           const numeric_index_type nnz,
                           const numeric_index_type noz,
                           const numeric_index_type blocksize_in)
{
  // So compilers don't warn when !LIBMESH_ENABLE_BLOCKED_STORAGE
  libmesh_ignore(blocksize_in);

  // Clear initialized matrices
  if (this->initialized())
    this->clear();

  this->_is_initialized = true;

  _vector = NumericVector<T>::build(communicator());

  _vector->init(m_in, m_l);
}


template <typename T>
void LumpedMatrix<T>::init (const numeric_index_type m_in,
                           const numeric_index_type n_in,
                           const numeric_index_type m_l,
                           const numeric_index_type n_l,
                           const std::vector<numeric_index_type> & n_nz,
                           const std::vector<numeric_index_type> & n_oz,
                           const numeric_index_type blocksize_in)
{
  libmesh_not_implemented();
}


template <typename T>
void LumpedMatrix<T>::init ()
{
  libmesh_assert(this->_dof_map);

  // Clear initialized matrices
  if (this->initialized())
    this->clear();

  this->_is_initialized = true;


  const numeric_index_type my_m = this->_dof_map->n_dofs();
  const numeric_index_type my_n = my_m;
  const numeric_index_type n_l  = this->_dof_map->n_dofs_on_processor(this->processor_id());
  const numeric_index_type m_l  = n_l;

  _vector = NumericVector<T>::build(communicator());

  _vector->init(my_m, m_l);
}


template <typename T>
void LumpedMatrix<T>::update_preallocation_and_zero ()
{
  libmesh_not_implemented();
}


template <typename T>
void LumpedMatrix<T>::zero ()
{
  libmesh_assert (this->initialized());

  _vector->zero();
}

template <typename T>
void LumpedMatrix<T>::zero_rows (std::vector<numeric_index_type> & rows, T diag_value)
{
  libmesh_assert (this->initialized());

  for (size_t i = 0; i < rows.size(); i++)
    _vector->set(rows[i]) = diag_value;
}

template <typename T>
void LumpedMatrix<T>::clear ()
{
  _vector->release();
}



template <typename T>
Real LumpedMatrix<T>::l1_norm () const
{
  return _vector->l1_norm();
}



template <typename T>
Real LumpedMatrix<T>::linfty_norm () const
{
  return _vector->linfty_norm();
}



template <typename T>
void LumpedMatrix<T>::print_matlab (const std::string & name) const
{
  libmesh_not_implemented();
}





template <typename T>
void LumpedMatrix<T>::print_personal(std::ostream & os) const
{
  libmesh_not_implemented();
}






template <typename T>
void LumpedMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                const std::vector<numeric_index_type> & rows,
                                const std::vector<numeric_index_type> & cols)
{
  const numeric_index_type n_rows = dm.m();
  const numeric_index_type n_cols = dm.n();

  libmesh_assert_equal_to (rows.size(), n_rows);
  libmesh_assert_equal_to (cols.size(), n_cols);

  for (auto i = 0; i < dm.m(); i++)
    for (auto j = 0; j < dm.n(); j++)
      _vector->add(rows[i], dm(i,j));
}






template <typename T>
void LumpedMatrix<T>::add_block_matrix(const DenseMatrix<T> & dm,
                                      const std::vector<numeric_index_type> & brows,
                                      const std::vector<numeric_index_type> & bcols)
{
  libmesh_not_implemented();
}





template <typename T>
void LumpedMatrix<T>::_get_submatrix(SparseMatrix<T> & submatrix,
                                    const std::vector<numeric_index_type> & rows,
                                    const std::vector<numeric_index_type> & cols,
                                    const bool reuse_submatrix) const
{
  libmesh_not_implemented();
}



template <typename T>
void LumpedMatrix<T>::get_diagonal (NumericVector<T> & dest) const
{
  dest = *_vector;
}



template <typename T>
void LumpedMatrix<T>::get_transpose (SparseMatrix<T> & dest) const
{
  libmesh_not_implemented();
}





template <typename T>
void LumpedMatrix<T>::close ()
{
  _vector->close();
}

template <typename T>
void LumpedMatrix<T>::flush ()
{
  _vector->flush();
}



template <typename T>
numeric_index_type LumpedMatrix<T>::m () const
{
  return _vector->size();
}



template <typename T>
numeric_index_type LumpedMatrix<T>::n () const
{
  return _n;
}



template <typename T>
numeric_index_type LumpedMatrix<T>::row_start () const
{
  return _vector->first_local_index();
}



template <typename T>
numeric_index_type LumpedMatrix<T>::row_stop () const
{
  return _vector->last_local_index();
}



template <typename T>
void LumpedMatrix<T>::set (const numeric_index_type i,
                          const numeric_index_type j,
                          const T value)
{
  libmesh_assert(i == j);

  _vector->set(i, value);
}



template <typename T>
void LumpedMatrix<T>::add (const numeric_index_type i,
                          const numeric_index_type j,
                          const T value)
{
  // This is how the lumping happens
  _vector->add(i, value);
}



template <typename T>
void LumpedMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                const std::vector<numeric_index_type> & dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}







template <typename T>
void LumpedMatrix<T>::add (const T a_in, SparseMatrix<T> & X_in)
{
  libmesh_not_implemented();
}




template <typename T>
T LumpedMatrix<T>::operator () (const numeric_index_type i_in,
                               const numeric_index_type j_in) const
{
  libmesh_assert(i_in == j_in);

  return (*_vector)(i_in);
}




template <typename T>
bool LumpedMatrix<T>::closed() const
{
  return _vector->closed();
}



template <typename T>
void LumpedMatrix<T>::swap(LumpedMatrix<T> & m_in)
{
  libmesh_not_implemented();
}



//------------------------------------------------------------------
// Explicit instantiations
template class LumpedMatrix<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
