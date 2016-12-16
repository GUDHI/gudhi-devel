/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016  INRIA
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
 * complex, that can be created from a `Rips_complex`.
 */
struct SimplicialComplexForRips {
  /** \brief Handle to specify the simplex filtration value. */
  typedef unspecified Filtration_value;

  /** \brief Inserts a given range `Gudhi::rips_complex::Rips_complex::OneSkeletonGraph` in the simplicial complex. */
  template<class OneSkeletonGraph>
  void insert_graph(const OneSkeletonGraph& skel_graph);

  /** \brief Expands the simplicial complex containing only its one skeleton until a given maximal dimension. */
  void expansion(int max_dim);

  /** \brief Returns the number of vertices in the simplicial complex. */
  std::size_t num_vertices();

};

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // CONCEPT_RIPS_COMPLEX_SIMPLICIAL_COMPLEX_FOR_RIPS_H_
