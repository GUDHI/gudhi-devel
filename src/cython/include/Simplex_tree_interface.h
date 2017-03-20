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
#include <gudhi/Points_off_io.h>

#include "Persistent_cohomology_interface.h"

#include <iostream>
#include <vector>
#include <utility>  // std::pair

namespace Gudhi {

template<typename SimplexTreeOptions = Simplex_tree_options_full_featured>
class Simplex_tree_interface : public Simplex_tree<SimplexTreeOptions> {
 public:
  using Base = Simplex_tree<SimplexTreeOptions>;
  using Filtration_value = typename Base::Filtration_value;
  using Vertex_handle = typename Base::Vertex_handle;
  using Simplex_handle = typename Base::Simplex_handle;
  using Insertion_result = typename std::pair<Simplex_handle, bool>;
  using Simplex = std::vector<Vertex_handle>;
  using Filtered_complex = std::pair<Simplex, Filtration_value>;
  using Complex_tree = std::vector<Filtered_complex>;

 public:

  bool find_simplex(const Simplex& vh) {
    return (Base::find(vh) != Base::null_simplex());
  }

  bool insert_simplex_and_subfaces(const Simplex& complex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex_and_subfaces(complex, filtration);
    Base::initialize_filtration();
    return (result.second);
  }

  // Do not interface this function, only used in strong witness interface for complex creation
  bool insert_simplex_and_subfaces(const std::vector<std::size_t>& complex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex_and_subfaces(complex, filtration);
    Base::initialize_filtration();
    return (result.second);
  }

  Filtration_value simplex_filtration(const Simplex& complex) {
    return Base::filtration(Base::find(complex));
  }

  void remove_maximal_simplex(const Simplex& complex) {
    Base::remove_maximal_simplex(Base::find(complex));
    Base::initialize_filtration();
  }

  Complex_tree get_filtered_tree() {
    Complex_tree filtered_tree;
    for (auto f_simplex : Base::filtration_simplex_range()) {
      Simplex simplex;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      filtered_tree.push_back(std::make_pair(simplex, Base::filtration(f_simplex)));
    }
    return filtered_tree;

  }

  Complex_tree get_skeleton_tree(int dimension) {
    Complex_tree skeleton_tree;
    for (auto f_simplex : Base::skeleton_simplex_range(dimension)) {
      Simplex simplex;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      skeleton_tree.push_back(std::make_pair(simplex, Base::filtration(f_simplex)));
    }
    return skeleton_tree;
  }

  Complex_tree get_star_tree(const Simplex& complex) {
    Complex_tree star_tree;
    for (auto f_simplex : Base::star_simplex_range(Base::find(complex))) {
      Simplex simplex;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      star_tree.push_back(std::make_pair(simplex, Base::filtration(f_simplex)));
    }
    return star_tree;
  }

  Complex_tree get_coface_tree(const Simplex& complex, int dimension) {
    Complex_tree coface_tree;
    for (auto f_simplex : Base::cofaces_simplex_range(Base::find(complex), dimension)) {
      Simplex simplex;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      coface_tree.push_back(std::make_pair(simplex, Base::filtration(f_simplex)));
    }
    return coface_tree;
  }

  void create_persistence(Gudhi::Persistent_cohomology_interface<Base>* pcoh) {
    pcoh = new Gudhi::Persistent_cohomology_interface<Base>(*this);
  }

};

} // namespace Gudhi

#endif  // SIMPLEX_TREE_INTERFACE_H

