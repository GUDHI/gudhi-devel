/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ZZ_STABLE_COMPLEX_TYPE_H_
#define CONCEPT_ZZ_STABLE_COMPLEX_TYPE_H_

/** @file StableFilteredComplex.h
 * @brief Contains @ref Gudhi::zigzag_persistence::StableFilteredComplex concept.
 */

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief Data structure storing the simplices and their filtration values in the current complex. 
 * The concept is realized for example by
 * @ref Gudhi::Simplex_tree < Gudhi::zigzag_persistence::Simplex_tree_options_oscillating_rips >.
 */
class StableFilteredComplex {
 public:
  /**
   * @brief Integer type that needs to be long enough to store the numbers of arrows in the zigzag filtration.
   */
  typename Simplex_key;

  /**
   * @brief Handle to specify a simplex. The simplex handles have to be stable, that is, they do not invalidate when
   * a simplex is added or removed from the complex (except for the removed simplices them-selves of course).
   */
  typename Simplex_handle;

  /**
   * @brief Type for filtration values. Usually 'double'.
   */
  typename Filtration_value;

  /**
   * @brief Removes the given simplex. Assumes that the simplex is maximal and can be safely removed.
   * 
   * @param sh Simplex handle representing the simplex to remove.
   */
  void remove_maximal_simplex(Simplex_handle sh);

  /**
   * @brief Adds a vertex or an edge in a flag complex, as well as all
   * simplices of its star, defined to maintain the property
   * of the complex to be a flag complex, truncated at dimension dim_max.
   *
   * @param u ID of one end of the edge.
   * @param v ID of the other end of the edge. If @p u == @p v, then the input is considered as a vertex.
   * @param fil Filtration value of the edge.
   * @param dim_max Maximal dimension of the expansion. If set to -1, the expansion goes as far as possible.
   * @param added_simplices Container for all new simplices induced by the insertion of the edge. 
   * If not empty at start, the content of the container should @b not be erased by the method.
   */
  void insert_edge_as_flag(int u,
                           int v,
                           Filtration_value fil,
                           int dim_max,
                           std::vector<Simplex_handle>& added_simplices);

  /**
   * @brief Returns the key associated to the given simplex.
   */
  Simplex_key key(Simplex_handle sh);

  /**
   * @brief Assignes the given value to the given simplex as a key.
   */
  void assign_key(Simplex_handle sh, Simplex_key key);

  /**
   * @brief Finds the given simplex in the complex and returns the associated simplex handle.
   * 
   * @tparam VertexRange Range over the vertices of a simplex.
   * @param simplex Simplex to find represented by its vertices.
   * @return The simplex handle associated to @p simplex if the simplex is found.
   */
  template <class VertexRange>
  Simplex_handle find(const VertexRange& simplex);

  /**
   * @brief Returns the filtration value of the given simplex.
   */
  Filtration_value filtration(Simplex_handle sh);

  /**
   * @brief Returns a range (with begin() and end() methods) over the star of the given simplex, including the simplex.
   * 
   * @param simplex Simplex to compute the star from.
   * @return A iterable range over the star of @p simplex.
   */
  auto star_simplex_range(const Simplex_handle simplex);

  /**
   * @brief Returns a range over the vertices of a simplex. The vertices have to be ordered monotonously by their
   * labels.
   */
  auto simplex_vertex_range(Simplex_handle sh) const;
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // CONCEPT_ZZ_STABLE_COMPLEX_TYPE_H_
