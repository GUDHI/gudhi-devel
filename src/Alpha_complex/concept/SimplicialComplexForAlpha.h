/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ALPHA_COMPLEX_SIMPLICIAL_COMPLEX_FOR_ALPHA_H_
#define CONCEPT_ALPHA_COMPLEX_SIMPLICIAL_COMPLEX_FOR_ALPHA_H_

namespace Gudhi {

namespace alpha_complex {

/** \brief The concept SimplicialComplexForAlpha describes the requirements for a type to implement a simplicial
 * complex, that can be created from a `Alpha_complex`.
 */
struct SimplicialComplexForAlpha {
  /** Handle to specify a simplex. */
  typedef unspecified Simplex_handle;
  /** Handle to specify a vertex. Must be a non-negative integer. */
  typedef unspecified Vertex_handle;
  /** Handle to specify the simplex filtration value. */
  typedef unspecified Filtration_value;

  /** Returns the number of vertices in the simplicial complex. */
  std::size_t num_vertices();

  /** Gets the 'simplex' dimension. */
  int dimension(Simplex_handle simplex);

  /** Assigns the 'simplex' with the given 'filtration' value. */
  int assign_filtration(Simplex_handle simplex, Filtration_value filtration);

  /** \brief Inserts a simplex with vertices from a given simplex (represented by a vector of Vertex_handle) in the
   * simplicial complex with the given 'filtration' value. */
  void insert_simplex_and_subfaces(std::vector<Vertex_handle> const & vertex_range, Filtration_value filtration);

  /** Browses the simplicial complex to make the filtration non-decreasing. */
  void make_filtration_non_decreasing();

  /** Prune the simplicial complex above 'filtration' value given as parameter. */
  void prune_above_filtration(Filtration_value filtration);

  /** \brief Iterator over vertices of a simplex.
   *
   * 'value type' must be 'Vertex_handle'.*/
  typedef unspecified Simplex_vertex_range;

  /** \brief Returns a range over vertices of a given simplex. */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle const & simplex);

  /** \brief Iterator over the boundaries of the complex, in an arbitrary order.
   *
   * 'value_type' must be 'Simplex_handle'.*/
  typedef unspecified Boundary_simplex_range;

  /** \brief Returns a range over boundaries of a given simplex. */
  Boundary_simplex_range boundary_simplex_range(Simplex_handle const & simplex);

  /** \brief Iterator over the simplices of the skeleton of the complex, for a given dimension.
   *
   * 'value_type' must be 'Simplex_handle'. */
  typedef unspecified Skeleton_simplex_range;
  /** \brief Returns a range over the simplices of the skeleton of the simplicial complex, for a given
   * dimension. */
  Skeleton_simplex_range skeleton_simplex_range;

  /** \brief Return type of an insertion of a simplex
   */
  typedef unspecified Insertion_result_type;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // CONCEPT_ALPHA_COMPLEX_SIMPLICIAL_COMPLEX_FOR_ALPHA_H_
