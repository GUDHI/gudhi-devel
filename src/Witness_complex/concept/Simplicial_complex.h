/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CONCEPT_WITNESS_COMPLEX_SIMPLICIAL_COMPLEX_H_
#define CONCEPT_WITNESS_COMPLEX_SIMPLICIAL_COMPLEX_H_

namespace Gudhi {

namespace witness_complex {

/** \brief The concept Simplicial_Complex describes the requirements 
 * for a type to implement a simplicial complex, 
 * used for example to build a 'Witness_complex'. 
 */
struct Simplicial_complex {
  /** Handle to specify a simplex. */
  typedef unspecified Simplex_handle;
  /** Handle to specify a vertex. Must be a non-negative integer. */
  typedef unspecified Vertex_handle;

  /** Returns a Simplex_hanlde that is different from all simplex handles 
   * of the simplices. */
  Simplex_handle null_simplex();

  /** \brief Iterator over the simplices of the complex,
   * in an arbitrary order.
   *
   * 'value_type' must be 'Simplex_handle'.*/
  typedef unspecified Complex_simplex_range;

  /**
   * \brief Returns a range over all the simplices of a
   * complex.
   */
  Complex_simplex_range complex_simplex_range();

  /** \brief Iterator over vertices of a simplex.
   *
   * 'value type' must be 'Vertex_handle'.*/
  typedef unspecified Simplex_vertex_range;

  /** \brief Returns a range over vertices of a given
   *  simplex. */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle const & simplex);

  /** \brief Return type of an insertion of a simplex
   */
  typedef unspecified Insertion_result_type;

  /** \brief Inserts a simplex with vertices from a given range
   *  'vertex_range' in the simplicial complex.
   *  */
  template< typedef Input_vertex_range >
  Insertion_result_type insert_simplex(Input_vertex_range const & vertex_range);

  /** \brief Finds a simplex with vertices given by a range
   *
   * If a simplex exists, its Simplex_handle is returned.
   * Otherwise null_simplex() is returned. */
  template< typedef Input_vertex_range >
  Simplex_handle find(Input_vertex_range const & vertex_range);
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // CONCEPT_WITNESS_COMPLEX_SIMPLICIAL_COMPLEX_H_
