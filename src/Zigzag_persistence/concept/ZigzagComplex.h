/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ZZ_COMPLEX_TYPE_H_
#define CONCEPT_ZZ_COMPLEX_TYPE_H_

/** @file ZigzagComplex.h
 * @brief Contains @ref Gudhi::zigzag_persistence::ZigzagComplex concept.
 */

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief Data structure storing the simplices in the current complex. 
 * The concept is realized for example by @ref Gudhi::Simplex_tree < Gudhi::Simplex_tree_options_zigzag_persistence >
 * or @ref Gudhi::Simplex_tree < Gudhi::Simplex_tree_options_zigzag_persistence_long >.
 */
class ZigzagComplex {
 public:
  /**
   * @brief Signed integer type that needs to be long enough to store the numbers of arrows in the zigzag filtration.
   */
  typename Simplex_key;

  /**
   * @brief Handle to specify a simplex.
   */
  typename Simplex_handle;

  /**
   * @brief Handle to specify a vertex. Should be an integer type.
   */
  typename Vertex_handle;

  /**
   * @brief Type for filtration values. Usually 'double'.
   */
  typename Filtration_value;

  /**
   * @brief Range of simplex handles over the boundary of a simplex
   */
  typename Boundary_simplex_range;

  /**
   * @brief Constructor
   */
  ZigzagComplex();

  /**
   * @brief Inserts the given simplex in the complex.
   * 
   * @tparam VertexRange Range over the vertices of a simplex.
   * @param simplex Simplex to insert represented by its vertices.
   * @return A pair of a simplex handle and a boolean. 
   * The simplex handle represents the inserted simplex and 
   * the boolean if simplex was already contained in the complex or not.
   */
  template <class VertexRange>
  std::pair<Simplex_handle, bool> insert_simplex(const VertexRange& simplex);

  /**
   * @brief Removes the given simplex. Assumes that the simplex is maximal and can be safely removed.
   * 
   * @param sh Simplex handle representing the simplex to remove.
   */
  void remove_maximal_simplex(Simplex_handle sh);

  /**
   * @brief Returns the dimension of the given simplex.
   * 
   * @param sh Simplex handle representing the simplex.
   * @return Dimension of @a sh.
   */
  int dimension(Simplex_handle sh);

  /**
   * @brief Returns the key associated to the given simplex.
   * 
   * @param sh Simplex handle representing the simplex.
   * @return The key.
   */
  Simplex_key key(Simplex_handle sh);

  /**
   * @brief Assignes the given value to the given simplex as a key.
   * 
   * @param sh Simplex handle representing the simplex.
   * @param key Values to associate as key.
   */
  void assign_key(Simplex_handle sh, Simplex_key key);

  /**
   * @brief Finds the given simplex in the complex and returns the associated simplex handle.
   * 
   * @tparam VertexRange Range over the vertices of a simplex.
   * @param simplex Simplex to find represented by its vertices.
   * @return The simplex handle associated to @a simplex if the simplex is found, @ref null_simplex() otherwise.
   */
  template <class VertexRange>
  Simplex_handle find(const VertexRange& simplex);

  /**
   * @brief Returns a range of simplex handles representing the boundary of the given simplex.
   * 
   * @param sh Simplex handle representing the simplex.
   * @return Range of simplex handles.
   */
  Boundary_simplex_range boundary_simplex_range(Simplex_handle sh);

  /**
   * @brief Returns a simplex handle representing a non existing simplex. 
   * 
   * @return A simplex handle.
   */
  Simplex_handle null_simplex();
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // CONCEPT_ZZ_COMPLEX_TYPE_H_
