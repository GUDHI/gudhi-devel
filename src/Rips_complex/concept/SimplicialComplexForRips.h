/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
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

#ifndef CONCEPT_RIPS_COMPLEX_SIMPLICIAL_COMPLEX_FOR_RIPS_H_
#define CONCEPT_RIPS_COMPLEX_SIMPLICIAL_COMPLEX_FOR_RIPS_H_

namespace Gudhi {

namespace rips_complex {

/** \brief The concept SimplicialComplexForRips describes the requirements for a type to implement a simplicial
 * complex, that can be created from a `Rips_complex`. The only available model for the moment is the `Simplex_tree`.
 */
struct SimplicialComplexForRips {
  /** \brief Type used to store the filtration values of the simplicial complex. */
  typedef unspecified Filtration_value;

  /** \brief Handle type to a simplex contained in the simplicial complex. */
  typedef unspecified Simplex_handle;

  /** \brief Inserts a given `Gudhi::rips_complex::Rips_complex::OneSkeletonGraph` in the simplicial complex. */
  template<class OneSkeletonGraph>
  void insert_graph(const OneSkeletonGraph& skel_graph);

  /** \brief Expands the simplicial complex containing only its one skeleton until a given maximal dimension as
   * explained in \ref ripsdefinition. */
  void expansion(int max_dim);

  /** \brief Expands a simplicial complex containing only a graph. Simplices corresponding to cliques in the graph are added
   * incrementally, faces before cofaces, unless the simplex has dimension larger than `max_dim` or `block_simplex`
   * returns true for this simplex.
   *
   * @param[in] max_dim Expansion maximal dimension value.
   * @param[in] block_simplex Blocker oracle. Its concept is <CODE>bool block_simplex(Simplex_handle sh)</CODE>
   *
   * The function identifies a candidate simplex whose faces are all already in the complex, inserts
   * it with a filtration value corresponding to the maximum of the filtration values of the faces, then calls
   * `block_simplex` on a `Simplex_handle` for this new simplex. If `block_simplex` returns true, the simplex is
   * removed, otherwise it is kept.
   */
  template< typename Blocker >
  void expansion_with_blockers(int max_dim, Blocker block_simplex);

  /** \brief Returns a range over the vertices of a simplex.  */
  unspecified simplex_vertex_range(Simplex_handle sh);

  /** \brief Returns the number of vertices in the simplicial complex. */
  std::size_t num_vertices();

};

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // CONCEPT_RIPS_COMPLEX_SIMPLICIAL_COMPLEX_FOR_RIPS_H_
