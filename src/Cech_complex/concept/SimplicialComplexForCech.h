/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
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

#ifndef CONCEPT_CECH_COMPLEX_SIMPLICIAL_COMPLEX_FOR_CECH_H_
#define CONCEPT_CECH_COMPLEX_SIMPLICIAL_COMPLEX_FOR_CECH_H_

namespace Gudhi {

namespace cech_complex {

/** \brief The concept SimplicialComplexForCech describes the requirements for a type to implement a simplicial
 * complex, that can be created from a `Cech_complex`.
 */
struct SimplicialComplexForCech {
  /** Handle to specify a simplex. */
  typedef unspecified Simplex_handle;
  /** Handle to specify a vertex. Must be a non-negative integer. */
  typedef unspecified Vertex_handle;
  /** Handle to specify the simplex filtration value. */
  typedef unspecified Filtration_value;

  /** Assigns the 'simplex' with the given 'filtration' value. */
  int assign_filtration(Simplex_handle simplex, Filtration_value filtration);

  /** \brief Returns a range over vertices of a given
   *  simplex. */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle const & simplex);

  /** \brief Inserts a given `Gudhi::ProximityGraph` in the simplicial complex. */
  template<class ProximityGraph>
  void insert_graph(const ProximityGraph& proximity_graph);

  /** \brief Expands the simplicial complex containing only its one skeleton until a given maximal dimension.
   * expansion can be blocked by the blocker oracle. */
  template< typename Blocker >
  void expansion_with_blockers(int max_dim, Blocker block_simplex);

  /** Returns the number of vertices in the simplicial complex. */
  std::size_t num_vertices();

};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // CONCEPT_ALPHA_COMPLEX_SIMPLICIAL_COMPLEX_FOR_ALPHA_H_
