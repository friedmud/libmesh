// $Id: face_tri6.h,v 1.16 2003-09-25 21:46:55 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __tri6_h__
#define __tri6_h__


// C++ includes

// Local includes
#include "libmesh_common.h"
#include "face_tri.h"



// Forward declarations




/**
 * The \p Tri6 is an element in 2D composed of 6 nodes.
 * It is numbered like this:
 * \verbatim
 *  TRI6:  2      
 *         o      
 *        / \     
 *       /   \    
 *    5 o     o 4 
 *     /       \  
 *    /         \ 
 *   o-----o-----o
 *   0     3     1
 * \endverbatim
 */

// ------------------------------------------------------------
// Tri6 class definition
class Tri6 : public Tri
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Tri6  (const Elem* p=NULL);

  /**
   * @returns \p TRI6
   */
  ElemType type ()   const { return TRI6; }

  /**
   * @returns 6
   */
  unsigned int n_nodes() const { return 6; }
  
  /**
   * @returns 4
   */
  unsigned int n_sub_elem() const { return 4; }
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   *
   * We reimplemenet this method here for the \p Quad8 since we can
   * use the center node of each edge to provide a perfect (unique)
   * key.
   */
  unsigned int key (const unsigned int s) const;
  
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sf=0) const;
  
    
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;

  unsigned int vtk_element_type (const unsigned int) const
  { return 6; }

  /**
   * @returns 2 for all \p n
   */
  unsigned int n_second_order_adjacent_vertices (const unsigned int) const
      { return 2; }

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   * Note that \p n is counted as depicted above, \f$ 3 \le n < 6 \f$.
   */
  unsigned short int second_order_adjacent_vertex (const unsigned int n,
						   const unsigned int v) const;
 
  
private:
  
  
#ifdef ENABLE_AMR
  
  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int i,
			  const unsigned int j,
			  const unsigned int k) const
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[4][6][6];
  
#endif


private:
  
  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[3][2];
    
};




// ------------------------------------------------------------
// Tri6 class member functions
inline
Tri6::Tri6(const Elem* p) :
  Tri(Tri6::n_nodes(), p) 
{
}


#endif
