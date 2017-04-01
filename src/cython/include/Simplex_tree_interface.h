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

#ifndef INCLUDE_SIMPLEX_TREE_INTERFACE_H_
#define INCLUDE_SIMPLEX_TREE_INTERFACE_H_

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
  using Complex = std::vector<std::pair<Simplex, Filtration_value>>;

 public:
  bool find_simplex(const Simplex& vh) {
    return (Base::find(vh) != Base::null_simplex());
  }

  bool insert(const Simplex& simplex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex_and_subfaces(simplex, filtration);
    return (result.second);
  }

  // Do not interface this function, only used in alpha complex interface for complex creation
  bool insert_simplex(const Simplex& simplex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex(simplex, filtration);
    return (result.second);
  }

  // Do not interface this function, only used in interface for complex creation
  bool insert_simplex_and_subfaces(const Simplex& simplex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex_and_subfaces(simplex, filtration);
    return (result.second);
  }

  // Do not interface this function, only used in strong witness interface for complex creation
  bool insert_simplex(const std::vector<std::size_t>& complex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex(complex, filtration);
    return (result.second);
  }

  // Do not interface this function, only used in strong witness interface for complex creation
  bool insert_simplex_and_subfaces(const std::vector<std::size_t>& complex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex_and_subfaces(complex, filtration);
    return (result.second);
  }

  Filtration_value simplex_filtration(const Simplex& simplex) {
    return Base::filtration(Base::find(simplex));
  }

  void remove_maximal_simplex(const Simplex& simplex) {
    Base::remove_maximal_simplex(Base::find(simplex));
    Base::initialize_filtration();
  }

  Complex get_filtration() {
    Base::initialize_filtration();
    Complex filtrations;
    for (auto f_simplex : Base::filtration_simplex_range()) {
      Simplex simplex;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      filtrations.push_back(std::make_pair(simplex, Base::filtration(f_simplex)));
    }
    return filtrations;
  }

  Complex get_skeleton(int dimension) {
    Complex skeletons;
    for (auto f_simplex : Base::skeleton_simplex_range(dimension)) {
      Simplex simplex;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      skeletons.push_back(std::make_pair(simplex, Base::filtration(f_simplex)));
    }
    return skeletons;
  }

  Complex get_star(const Simplex& simplex) {
    Complex star;
    for (auto f_simplex : Base::star_simplex_range(Base::find(simplex))) {
      Simplex simplex_star;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        std::cout << vertex << " ";
        simplex_star.insert(simplex_star.begin(), vertex);
      }
      std::cout << std::endl;
      star.push_back(std::make_pair(simplex_star, Base::filtration(f_simplex)));
    }
    return star;
  }

  Complex get_cofaces(const Simplex& simplex, int dimension) {
    Complex cofaces;
    for (auto f_simplex : Base::cofaces_simplex_range(Base::find(simplex), dimension)) {
      Simplex simplex_coface;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        std::cout << vertex << " ";
        simplex_coface.insert(simplex_coface.begin(), vertex);
      }
      std::cout << std::endl;
      cofaces.push_back(std::make_pair(simplex_coface, Base::filtration(f_simplex)));
    }
    return cofaces;
  }

  void create_persistence(Gudhi::Persistent_cohomology_interface<Base>* pcoh) {
    Base::initialize_filtration();
    pcoh = new Gudhi::Persistent_cohomology_interface<Base>(*this);
  }
};

}  // namespace Gudhi

#endif  // INCLUDE_SIMPLEX_TREE_INTERFACE_H_
