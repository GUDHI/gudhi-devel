/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_SIMPLEX_TREE_OPTIONS_H_
#define SIMPLEX_TREE_SIMPLEX_TREE_OPTIONS_H_

#include <gudhi/Simplex_tree/indexing_tag.h>

#include <cstdint>

namespace Gudhi {

/** \addtogroup simplex_tree 
 * Pre-defined options for the Simplex_tree.
 *  @{
 */

/** Model of SimplexTreeOptions.
 * 
 * Default options version of the Simplex_tree.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */
struct Simplex_tree_options_default {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  static const bool link_nodes_by_label = false;
  static const bool stable_simplex_handles = false;
  static const bool is_multi_parameter = false;
};

/** Model of SimplexTreeOptions.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */
struct Simplex_tree_options_full_featured {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  static const bool link_nodes_by_label = true;
  static const bool stable_simplex_handles = true;
  static const bool is_multi_parameter = false;
};

/** Model of SimplexTreeOptions.
 * 
 * Minimal version of the Simplex_tree. No filtration values are stored and it is impossible to compute persistence
 * with these options. */
struct Simplex_tree_options_minimal {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = false;
  static const bool store_filtration = false;
  static const bool contiguous_vertices = false;
  static const bool link_nodes_by_label = false;
  static const bool stable_simplex_handles = false;
  static const bool is_multi_parameter = false;
};

/** @private @brief Model of SimplexTreeOptions, faster than `Simplex_tree_options_default` but note the unsafe
 * `contiguous_vertices` option.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */
struct Simplex_tree_options_fast_persistence {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef float Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = true;
  static const bool link_nodes_by_label = false;
  static const bool stable_simplex_handles = false;
  static const bool is_multi_parameter = false;
};

/** @}*/  // end addtogroup simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SIMPLEX_TREE_OPTIONS_H_
