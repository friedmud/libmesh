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



#ifndef LIBMESH_LUMPED_MATRIX_H
#define LIBMESH_LUMPED_MATRIX_H

#include "libmesh/libmesh_common.h"

// Local includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/parallel.h"
#include "libmesh/numeric_vector.h"

// C++ includes
#include <algorithm>

// Macro to identify and debug functions which should be called in
// parallel on parallel matrices but which may be called in serial on
// serial matrices.  This macro will only be valid inside non-static
// LumpedMatrix methods
#undef semiparallel_only
#ifndef NDEBUG
#include <cstring>

#define semiparallel_only() do { if (this->initialized()) { const char * mytype; \
      MatGetType(_mat,&mytype);                                         \
      if (!strcmp(mytype, MATSEQAIJ))                                   \
        parallel_object_only(); } } while (0)
#else
#define semiparallel_only()
#endif


namespace libMesh
{

// Forward Declarations
template <typename T> class DenseMatrix;


/**
 * This class provides a SparseMatrix looking object that
 * actually only stores a diagonal and automatically does "lumping"
 * (every off-diagonal entry is just added to the diagonal)
 *
 * It is implemented as a NumericVector underneath
 *
 * \author Derek Gaston
 * \date 2018
 * \brief Diagonal lumped matrix
 */
template <typename T>
class LumpedMatrix libmesh_final : public SparseMatrix<T>
{
public:
  /**
   * Constructor; initializes the matrix to be empty, without any
   * structure, i.e.  the matrix is not usable at all. This
   * constructor is therefore only useful for matrices which are
   * members of a class. All other matrices should be created at a
   * point in the data flow where all necessary information is
   * available.
   *
   * You have to initialize the matrix before usage with \p init(...).
   */
  explicit
  LumpedMatrix (const Parallel::Communicator & comm_in
               LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor. Free all memory, but do not release the memory of the
   * sparsity structure.
   */
  ~LumpedMatrix ();

  /**
   * Note, only m,n and m_l are actually used.
   */
  virtual void init (const numeric_index_type m,
                     const numeric_index_type n,
                     const numeric_index_type m_l,
                     const numeric_index_type n_l,
                     const numeric_index_type nnz=30,
                     const numeric_index_type noz=10,
                     const numeric_index_type blocksize=1) libmesh_override;

  /**
   * Note, only m,n and m_l are actually used.
   *
   * The rest of the interface is kept for consistency.
   *
   * \param m The global number of rows.
   * \param n The global number of columns.
   * \param m_l The local number of rows.
   * \param n_l The local number of columns.
   * \param n_nz array containing the number of nonzeros in each row of the DIAGONAL portion of the local submatrix.
   * \param n_oz Array containing the number of nonzeros in each row of the OFF-DIAGONAL portion of the local submatrix.
   * \param blocksize Optional value indicating dense coupled blocks for systems with multiple variables all of the same type.
   */
  void init (const numeric_index_type m,
             const numeric_index_type n,
             const numeric_index_type m_l,
             const numeric_index_type n_l,
             const std::vector<numeric_index_type> & n_nz,
             const std::vector<numeric_index_type> & n_oz,
             const numeric_index_type blocksize=1);

  virtual void init () libmesh_override;

  /**
   * Update the sparsity pattern based on \p dof_map, and set the matrix
   * to zero. This is useful in cases where the sparsity pattern changes
   * during a computation.
   */
  void update_preallocation_and_zero();

  virtual void clear () libmesh_override;

  virtual void zero () libmesh_override;

  virtual void zero_rows (std::vector<numeric_index_type> & rows, T diag_value = 0.0) libmesh_override;

  virtual void close () libmesh_override;

  virtual void flush () libmesh_override;

  virtual numeric_index_type m () const libmesh_override;

  virtual numeric_index_type n () const libmesh_override;

  virtual numeric_index_type row_start () const libmesh_override;

  virtual numeric_index_type row_stop () const libmesh_override;

  virtual void set (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) libmesh_override;

  virtual void add (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) libmesh_override;

  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & rows,
                           const std::vector<numeric_index_type> & cols) libmesh_override;

  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & dof_indices) libmesh_override;

  virtual void add_block_matrix (const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & brows,
                                 const std::vector<numeric_index_type> & bcols) libmesh_override;

  virtual void add_block_matrix (const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & dof_indices) libmesh_override
  { this->add_block_matrix (dm, dof_indices, dof_indices); }

  /**
   * Compute A += a*X for scalar \p a, matrix \p X.
   *
   * \note The matrices A and X need to have the same nonzero pattern,
   *
   * \note It is advisable to not only allocate appropriate memory
   * with \p init(), but also explicitly zero the terms of \p this
   * whenever you add a non-zero value to \p X.
   *
   * \note \p X will be closed, if not already done, before performing
   * any work.
   */
  virtual void add (const T a, SparseMatrix<T> & X) libmesh_override;

  virtual T operator () (const numeric_index_type i,
                         const numeric_index_type j) const libmesh_override;

  virtual Real l1_norm () const libmesh_override;

  virtual Real linfty_norm () const libmesh_override;

  virtual bool closed() const libmesh_override;

  /**
   * Print the contents of the matrix to the screen with the PETSc
   * viewer. This function only allows printing to standard out since
   * we have limited ourselves to one PETSc implementation for
   * writing.
   */
  virtual void print_personal(std::ostream & os=libMesh::out) const libmesh_override;

  virtual void print_matlab(const std::string & name = "") const libmesh_override;

  virtual void get_diagonal (NumericVector<T> & dest) const libmesh_override;

  virtual void get_transpose (SparseMatrix<T> & dest) const libmesh_override;

  /**
   * Swaps the internal data pointers of two PetscMatrices, no actual
   * values are swapped.
   */
  void swap (LumpedMatrix<T> &);

protected:
  std::unique_ptr<NumericVector<T>> _vector;
};

} // namespace libMesh

#endif
