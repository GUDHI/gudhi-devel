/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_CECH_COMPLEX_SIMPLICIAL_COMPLEX_FOR_MEB_H_
#define CONCEPT_CECH_COMPLEX_SIMPLICIAL_COMPLEX_FOR_MEB_H_

namespace Gudhi {

namespace cech_complex {

/** The concept SimplicialComplexForMEB describes the requirements for a type to implement a simplicial
 * complex that can be filled by `assign_MEB_filtration()`. It is typically satisfied by `Simplex_tree` if
 * `SimplexTreeOptions::store_key` is true.
 */
struct SimplicialComplexForMEB {
  /** \brief Handle for a simplex. */
  typedef unspecified Simplex_handle;
  /** \brief Handle for a vertex. Must be a non-negative integer,
   * it is also used as an index into the input list of points. */
  typedef unspecified Vertex_handle;
  /** \brief Type of filtration values. */
  typedef unspecified Filtration_value;
  /** \brief Integer type large enough to index all simplices. */
  typedef unspecified Simplex_key;

  /** \brief Returns the filtration value to the 'simplex'. */
  Filtration_value filtration(Simplex_handle simplex);
  /** \brief Assigns this 'filtration' value to the 'simplex'. */
  int assign_filtration(Simplex_handle simplex, Filtration_value filtration);

  /** \brief Returns the key assigned to the 'simplex' with `assign_key()`. */
  Simplex_key key(Simplex_handle simplex);
  /** \brief Assigns this 'key' to the 'simplex'. */
  void assign_key(Simplex_handle simplex, Simplex_key key);

  /** \brief Returns a range over vertices (as Vertex_handle) of a given simplex. */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle simplex);

  /** \brief Returns a range of the pairs (simplex, opposite vertex) of the boundary of the 'simplex'. */
  Boundary_opposite_vertex_simplex_range boundary_opposite_vertex_simplex_range(Simplex_handle simplex);

  /** \brief Calls `callback(simplex, dim)` for every simplex of the complex,
   * with the guarantee that faces are visited before cofaces. */
  void for_each_simplex(auto callback);
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CONCEPT_CECH_COMPLEX_SIMPLICIAL_COMPLEX_FOR_MEB_H_
