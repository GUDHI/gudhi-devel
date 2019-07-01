/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ALPHA_COMPLEX_SIMPLICIAL_COMPLEX_FOR_ALPHA_3D_H_
#define CONCEPT_ALPHA_COMPLEX_SIMPLICIAL_COMPLEX_FOR_ALPHA_3D_H_

namespace Gudhi {

namespace alpha_complex {

/** \brief The concept SimplicialComplexForAlpha3d describes the requirements for a type to implement a simplicial
 * complex, that can be created from a `Alpha_complex_3d`.
 */
struct SimplicialComplexForAlpha3d {
  /** Handle to specify a vertex. Must be a non-negative integer. */
  typedef unspecified Vertex_handle;
  /** Handle to specify the simplex filtration value. */
  typedef unspecified Filtration_value;

  /** Returns the number of vertices in the simplicial complex. */
  std::size_t num_vertices();

  /** \brief Inserts a simplex from a given simplex (represented by a vector of Vertex_handle) in the
   * simplicial complex with the given 'filtration' value. */
  void insert_simplex(std::vector<Vertex_handle> const& vertex_range, Filtration_value filtration);

  /** Browses the simplicial complex to make the filtration non-decreasing. */
  void make_filtration_non_decreasing();

  /** Prune the simplicial complex above 'filtration' value given as parameter. */
  void prune_above_filtration(Filtration_value filtration);
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // CONCEPT_ALPHA_COMPLEX_SIMPLICIAL_COMPLEX_FOR_ALPHA_3D_H_
