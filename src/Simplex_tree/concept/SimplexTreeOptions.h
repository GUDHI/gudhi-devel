 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       Marc Glisse
  *
  *    Copyright (C) 2015  INRIA Saclay - Ile-de-France (France)
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
  /// Must be a signed integer type.
  typedef SimplexKey Simplex_key;
  /// If true, each simplex has extra storage for one `Simplex_key`. Necessary for `Persistent_cohomology`.
  static const bool store_key;
  /// If true, each simplex has extra storage for one `Filtration_value`, and this value is propagated by operations like `Gudhi::Simplex_tree::expansion`. Without it, `Persistent_cohomology` degenerates to computing usual (non-persistent) cohomology.
  static const bool store_filtration;
  /// If true, the list of vertices present in the complex must always be 0, ..., num_vertices-1, without any hole.
  static constexpr bool contiguous_vertices;
};

