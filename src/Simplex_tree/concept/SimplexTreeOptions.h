/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
  
/** \brief Concept of the template parameter for the class `Gudhi::Simplex_tree<SimplexTreeOptions>`.
 *
 * One model for this is `Gudhi::Simplex_tree_options_full_featured`. If you want to provide your own, it is recommended that you derive from it and override some parts instead of writing a class from scratch.
 */
struct SimplexTreeOptions {
  /// Forced for now.
  typedef IndexingTag Indexing_tag;
  /// Must be a signed integer type. It admits a total order <.
  typedef VertexHandle Vertex_handle;
  /// Must be comparable with operator<.
  typedef FiltrationValue Filtration_value;
  /// Must be an integer type.
  typedef SimplexKey Simplex_key;
  /// If true, each simplex has extra storage for one `Simplex_key`. Necessary for `Persistent_cohomology`.
  static const bool store_key;
  /// If true, each simplex has extra storage for one `Filtration_value`, and this value is propagated by operations like `Gudhi::Simplex_tree::expansion`. Without it, `Persistent_cohomology` degenerates to computing usual (non-persistent) cohomology.
  static const bool store_filtration;
  /// If true, the list of vertices present in the complex must always be 0, ..., num_vertices-1, without any hole.
  static constexpr bool contiguous_vertices;
  /// If true, the lists of `Node` with same label are stored to enhance cofaces and stars access.
  static const bool link_nodes_by_label;
  /// If true, Simplex_handle will not be invalidated after insertions or removals.
  static const bool stable_simplex_handles;
  /// If true, assumes that Filtration_value is vector-like instead of float-like. This also assumes that Filtration_values is a class, which has a push_to method to push a filtration value $x$ onto $this>=0$. 
  static const bool is_multi_parameter;
};

