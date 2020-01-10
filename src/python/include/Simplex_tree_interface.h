/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_SIMPLEX_TREE_INTERFACE_H_
#define INCLUDE_SIMPLEX_TREE_INTERFACE_H_

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/Simplex_tree/Simplex_tree_iterators.h>

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
  using Filtered_simplices = std::vector<std::pair<Simplex, Filtration_value>>;

 public:
  bool find_simplex(const Simplex& vh) {
    return (Base::find(vh) != Base::null_simplex());
  }

  void assign_simplex_filtration(const Simplex& vh, Filtration_value filtration) {
    Base::assign_filtration(Base::find(vh), filtration);
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
  bool insert_simplex(const std::vector<std::size_t>& simplex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex(simplex, filtration);
    return (result.second);
  }

  // Do not interface this function, only used in strong witness interface for complex creation
  bool insert_simplex_and_subfaces(const std::vector<std::size_t>& simplex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex_and_subfaces(simplex, filtration);
    return (result.second);
  }

  Filtration_value simplex_filtration(const Simplex& simplex) {
    return Base::filtration(Base::find(simplex));
  }

  void remove_maximal_simplex(const Simplex& simplex) {
    Base::remove_maximal_simplex(Base::find(simplex));
    Base::initialize_filtration();
  }

  // If called with key < 0: allocate a new key for current iteration.
  // If called with key >= 0: remove the key from used set (iteration is over)
  int init_or_stop_iteration(int key) {
    static bool filtration_initialized = false;
    if (key < 0 && !filtration_initialized)
    {
      Base::initialize_filtration();
      filtration_initialized = true;
    }
    static std::set<int> used_keys;
    if (key >= 0)
      used_keys.erase(key);
    else
    {
      // Search for smallest int available and use it:
      key = 0;
      while (used_keys.count(key))
        key++;
    }
    return key;
  }

  std::pair<std::vector<int>, double> get_next_in_filtration(int key) {
    static std::map<int, Simplex_tree_complex_simplex_iterator<Base>> iters;
    if (!iters.count(key))
      iters[key] = Simplex_tree_complex_simplex_iterator<Base>(this);
    if (iters[key] == Simplex_tree_complex_simplex_iterator<Base>())
    {
      iters.erase(key); //iteration is over
      init_or_stop_iteration(key); //"stop" because key >= 0
      return std::make_pair<Simplex, double>(std::vector<int>(), 0.0);
    }
    auto res = *(iters[key]++);
    Simplex simplex;
    for (auto vertex : Base::simplex_vertex_range(res))
      simplex.insert(simplex.begin(), vertex);
    return std::make_pair<Simplex, double>(std::move(simplex), res.get_ptr()->first);
  }

  Filtered_simplices get_filtration() {
    Base::initialize_filtration();
    Filtered_simplices filtrations;
    for (auto f_simplex : Base::filtration_simplex_range()) {
      Simplex simplex;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      filtrations.push_back(std::make_pair(simplex, Base::filtration(f_simplex)));
    }
    return filtrations;
  }

  Filtered_simplices get_skeleton(int dimension) {
    Filtered_simplices skeletons;
    for (auto f_simplex : Base::skeleton_simplex_range(dimension)) {
      Simplex simplex;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      skeletons.push_back(std::make_pair(simplex, Base::filtration(f_simplex)));
    }
    return skeletons;
  }

  Filtered_simplices get_star(const Simplex& simplex) {
    Filtered_simplices star;
    for (auto f_simplex : Base::star_simplex_range(Base::find(simplex))) {
      Simplex simplex_star;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex_star.insert(simplex_star.begin(), vertex);
      }
      star.push_back(std::make_pair(simplex_star, Base::filtration(f_simplex)));
    }
    return star;
  }

  Filtered_simplices get_cofaces(const Simplex& simplex, int dimension) {
    Filtered_simplices cofaces;
    for (auto f_simplex : Base::cofaces_simplex_range(Base::find(simplex), dimension)) {
      Simplex simplex_coface;
      for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
        simplex_coface.insert(simplex_coface.begin(), vertex);
      }
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
