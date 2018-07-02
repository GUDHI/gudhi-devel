/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2014 Inria
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

#ifndef CONCEPT_WITNESS_COMPLEX_SIMPLICIAL_COMPLEX_FOR_WITNESS_H_
#define CONCEPT_WITNESS_COMPLEX_SIMPLICIAL_COMPLEX_FOR_WITNESS_H_

namespace Gudhi {

namespace witness_complex {

/** \brief The concept SimplicialComplexForWitness describes the requirements 
 * for a type to implement a simplicial complex, 
 * used for example to build a Witness_complex or Strong_witness_complex. 
 */
struct SimplicialComplexForWitness {
  /** Handle to specify a simplex. */
  typedef unspecified Simplex_handle;
  // /** Handle to specify a vertex. Must be a non-negative integer. */
  // typedef unspecified Vertex_handle;

  /** \brief Returns a Simplex_hanlde that is different from all simplex handles 
   * of the simplices. */
  Simplex_handle null_simplex();

  /** Returns the number of vertices in the simplicial complex
   */
  std::size_t num_vertices();
  
  /** \brief Return type of an insertion of a simplex
   */
  typedef unspecified Insertion_result_type;

  /** \brief Inserts a simplex with vertices from a given range
   *  'vertex_range' in the simplicial complex.
   *  The function is only used in Witness_complex class
   *  and by construction, it is not necessary to check if 
   *  the faces are in the simplicial complex before insertion. 
   *  The simplex is given the filtration value 'filtration'.
   *  Filtration_value should be convertible from double.
   *  The return type is not used.
   *  */
  template< typedef Input_vertex_range >
  Insertion_result_type insert_simplex(Input_vertex_range const & vertex_range, Filtration_value filtration);

  /** \brief Inserts a simplex and all its faces
   *  with vertices from a given range
   *  'vertex_range' in the simplicial complex.
   *  The function is only used in Strong_witness_complex class.
   *  All inserted simplices are given the filtration
   *  value 'filtration'.
   *  Filtration_value should be convertible from double.
   *  The return type is not used.
   */

  template< typedef Input_vertex_range,
            typedef Filtration_value>
  Insertion_result_type insert_simplex_and_subfaces(Input_vertex_range const & vertex_range, Filtration_value filtration);
  
  /** \brief Finds a simplex with vertices given by a range
   *
   * If a simplex exists, its Simplex_handle is returned.
   * Otherwise null_simplex() is returned. */
  template< typedef Input_vertex_range >
  Simplex_handle find(Input_vertex_range const & vertex_range);

  /** \brief Sets the dimension of the simplicial complex to 
   * 'dimension'.
   */
  void set_dimension(int dimension);

  /** \brief Returns the filtration of the simplex given by
   *  the simplex handle 'sh'.
   */
  double filtration(Simplex_handle sh);
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // CONCEPT_WITNESS_COMPLEX_SIMPLICIAL_COMPLEX_FOR_WITNESS_H_
