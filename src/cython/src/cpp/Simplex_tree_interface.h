/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef SIMPLEX_TREE_INTERFACE_H
#define	SIMPLEX_TREE_INTERFACE_H

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>

#include <vector>
#include <utility>  // std::pair
#include <iostream>

namespace Gudhi {

template<typename SimplexTreeOptions = Simplex_tree_options_full_featured>
class Simplex_tree_interface : public Simplex_tree<SimplexTreeOptions> {
  typedef typename Simplex_tree<SimplexTreeOptions>::Simplex_handle Simplex_handle;
  typedef typename std::pair<Simplex_handle, bool> Insertion_result;
  using Simplex = std::vector<Vertex_handle>;
  using Filtered_complex = std::pair<Simplex, Filtration_value>;
  using Complex_tree = std::vector<Filtered_complex>;

 public:

  bool find_simplex(const Simplex& vh) {
    return (Simplex_tree<SimplexTreeOptions>::find(vh) != Simplex_tree<SimplexTreeOptions>::null_simplex());
  }

  bool insert_simplex_and_subfaces(const Simplex& complex, Filtration_value filtration = 0) {
    Insertion_result result = Simplex_tree<SimplexTreeOptions>::insert_simplex_and_subfaces(complex, filtration);
    return (result.second);
  }

  Filtration_value simplex_filtration(const Simplex& complex) {
    return Simplex_tree<SimplexTreeOptions>::filtration(Simplex_tree<SimplexTreeOptions>::find(complex));
  }

  void remove_maximal_simplex(const Simplex& complex) {
    return Simplex_tree<SimplexTreeOptions>::remove_maximal_simplex(Simplex_tree<SimplexTreeOptions>::find(complex));
  }

  Complex_tree get_filtered_tree() {
    Complex_tree filtered_tree;
    for (auto f_simplex : Simplex_tree<SimplexTreeOptions>::filtration_simplex_range()) {
      Simplex simplex;
      for (auto vertex : Simplex_tree<SimplexTreeOptions>::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      filtered_tree.push_back(std::make_pair(simplex, Simplex_tree<SimplexTreeOptions>::filtration(f_simplex)));
    }
    return filtered_tree;

  }

  Complex_tree get_skeleton_tree(int dimension) {
    Complex_tree skeleton_tree;
    for (auto f_simplex : Simplex_tree<SimplexTreeOptions>::skeleton_simplex_range(dimension)) {
      Simplex simplex;
      for (auto vertex : Simplex_tree<SimplexTreeOptions>::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      skeleton_tree.push_back(std::make_pair(simplex, Simplex_tree<SimplexTreeOptions>::filtration(f_simplex)));
    }
    return skeleton_tree;
  }

  Complex_tree get_star_tree(const Simplex& complex) {
    Complex_tree star_tree;
    for (auto f_simplex : Simplex_tree<SimplexTreeOptions>::star_simplex_range(Simplex_tree<SimplexTreeOptions>::find(complex))) {
      Simplex simplex;
      for (auto vertex : Simplex_tree<SimplexTreeOptions>::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      star_tree.push_back(std::make_pair(simplex, Simplex_tree<SimplexTreeOptions>::filtration(f_simplex)));
    }
    return star_tree;
  }

  Complex_tree get_coface_tree(const Simplex& complex, int dimension) {
    Complex_tree coface_tree;
    for (auto f_simplex : Simplex_tree<SimplexTreeOptions>::cofaces_simplex_range(Simplex_tree<SimplexTreeOptions>::find(complex), dimension)) {
      Simplex simplex;
      for (auto vertex : Simplex_tree<SimplexTreeOptions>::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      coface_tree.push_back(std::make_pair(simplex, Simplex_tree<SimplexTreeOptions>::filtration(f_simplex)));
    }
    return coface_tree;
  }

  void graph_expansion(std::vector<std::vector<double>>&points, int max_dimension, double max_edge_length) {
    Graph_t prox_graph = compute_proximity_graph(points, max_edge_length, euclidean_distance<std::vector<double>>);
    Simplex_tree<SimplexTreeOptions>::insert_graph(prox_graph);
    Simplex_tree<SimplexTreeOptions>::expansion(max_dimension);
    Simplex_tree<SimplexTreeOptions>::initialize_filtration();
  }

};

struct Simplex_tree_options_mini : Simplex_tree_options_full_featured {
  // Not doing persistence, so we don't need those
  static const bool store_key = true;
  static const bool store_filtration = false;
  // I have few vertices
  typedef short Vertex_handle;
};

} // namespace Gudhi

#endif  // SIMPLEX_TREE_INTERFACE_H

