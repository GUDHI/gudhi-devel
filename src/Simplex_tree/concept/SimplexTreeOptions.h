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
  /// Must be an integer type. Must be signed if store_morse_pairing == true.
  typedef SimplexKey Simplex_key;
  /// If true, each simplex has extra storage for one `Simplex_key`. Necessary for `Persistent_cohomology`.
  static const bool store_key;
  /// If true, each simplex has extra storage for one `Filtration_value`, and this value is propagated by operations like `Gudhi::Simplex_tree::expansion`. Without it, `Persistent_cohomology` degenerates to computing usual (non-persistent) cohomology.
  static const bool store_filtration;
  /// If true, the list of vertices present in the complex must always be 0, ..., num_vertices-1, without any hole.
  static constexpr bool contiguous_vertices;
  /// If true, allows insertion and deletion of simplices without invalidating iterators (Simplex_handles).
  static const bool simplex_handle_strong_validity = true;
  /// If true, stores two extra pointers in each Node. All Nodes with a same label u are linked this way in an intrusive list. 
  static const bool link_nodes_by_label = true;
  /// If true, stores an extra Simplex_handle in each Node to encode a pairing between simplices, such as a Morse matching. Simplex_key must be signed in that case.
  static const bool store_morse_matching = true;
};

