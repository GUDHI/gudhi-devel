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
#include <gudhi/Flag_complex_edge_collapser.h>

#include <iostream>
#include <vector>
#include <utility>  // std::pair
#include <tuple>
#include <iterator>  // for std::distance

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
  using Simplex_and_filtration = std::pair<Simplex, Filtration_value>;
  using Filtered_simplices = std::vector<Simplex_and_filtration>;
  using Skeleton_simplex_iterator = typename Base::Skeleton_simplex_iterator;
  using Complex_simplex_iterator = typename Base::Complex_simplex_iterator;
  using Extended_filtration_data = typename Base::Extended_filtration_data;
  typedef bool (*blocker_func)(Simplex simplex, void *user_data);

 public:

  Extended_filtration_data efd;
  
  bool find_simplex(const Simplex& simplex) {
    return (Base::find(simplex) != Base::null_simplex());
  }

  void assign_simplex_filtration(const Simplex& simplex, Filtration_value filtration) {
    Base::assign_filtration(Base::find(simplex), filtration);
    Base::clear_filtration();
  }

  bool insert(const Simplex& simplex, Filtration_value filtration = 0) {
    Insertion_result result = Base::insert_simplex_and_subfaces(simplex, filtration);
    if (result.first != Base::null_simplex())
      Base::clear_filtration();
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
    Base::clear_filtration();
  }

  Simplex_and_filtration get_simplex_and_filtration(Simplex_handle f_simplex) {
    Simplex simplex;
    for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
      simplex.insert(simplex.begin(), vertex);
    }
    return std::make_pair(std::move(simplex), Base::filtration(f_simplex));
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

  void compute_extended_filtration() {
    this->efd = this->extend_filtration();
    return;
  }

  std::vector<std::vector<std::pair<int, std::pair<Filtration_value, Filtration_value>>>> compute_extended_persistence_subdiagrams(const std::vector<std::pair<int, std::pair<Filtration_value, Filtration_value>>>& dgm, Filtration_value min_persistence){
    std::vector<std::vector<std::pair<int, std::pair<Filtration_value, Filtration_value>>>> new_dgm(4);
    for (unsigned int i = 0; i < dgm.size(); i++){
      std::pair<Filtration_value, Extended_simplex_type> px = this->decode_extended_filtration(dgm[i].second.first, this->efd);
      std::pair<Filtration_value, Extended_simplex_type> py = this->decode_extended_filtration(dgm[i].second.second, this->efd);
      std::pair<int, std::pair<Filtration_value, Filtration_value>> pd_point = std::make_pair(dgm[i].first, std::make_pair(px.first, py.first));
      if(std::abs(px.first - py.first) > min_persistence){
        //Ordinary
        if (px.second == Extended_simplex_type::UP && py.second == Extended_simplex_type::UP){
          new_dgm[0].push_back(pd_point);
        }
        // Relative
        else if (px.second == Extended_simplex_type::DOWN && py.second == Extended_simplex_type::DOWN){
          new_dgm[1].push_back(pd_point);
        }
        else{
          // Extended+
          if (px.first < py.first){
            new_dgm[2].push_back(pd_point);
          }
          //Extended-
          else{
            new_dgm[3].push_back(pd_point);
          }
        }
      }
    }
    return new_dgm;
  }

  Simplex_tree_interface* collapse_edges(int nb_collapse_iteration) {
    using Filtered_edge = std::tuple<Vertex_handle, Vertex_handle, Filtration_value>;
    std::vector<Filtered_edge> edges;
    for (Simplex_handle sh : Base::skeleton_simplex_range(1)) {
      if (Base::dimension(sh) == 1) {
        typename Base::Simplex_vertex_range rg = Base::simplex_vertex_range(sh);
        auto vit = rg.begin();
        Vertex_handle v = *vit;
        Vertex_handle w = *++vit;
        edges.emplace_back(v, w, Base::filtration(sh));
      }
    }

    for (int iteration = 0; iteration < nb_collapse_iteration; iteration++) {
      edges = Gudhi::collapse::flag_complex_collapse_edges(edges);
    }
    Simplex_tree_interface* collapsed_stree_ptr = new Simplex_tree_interface();
    // Copy the original 0-skeleton
    for (Simplex_handle sh : Base::skeleton_simplex_range(0)) {
      collapsed_stree_ptr->insert({*(Base::simplex_vertex_range(sh).begin())}, Base::filtration(sh));
    }
    // Insert remaining edges
    for (auto remaining_edge : edges) {
      collapsed_stree_ptr->insert({std::get<0>(remaining_edge), std::get<1>(remaining_edge)}, std::get<2>(remaining_edge));
    }
    return collapsed_stree_ptr;
  }

  void expansion_with_blockers_callback(int dimension, blocker_func user_func, void *user_data) {
    Base::expansion_with_blockers(dimension, [&](Simplex_handle sh){
      Simplex simplex;
      for (auto vertex : Base::simplex_vertex_range(sh)) {
        simplex.insert(simplex.begin(), vertex);
      }
      return user_func(simplex, user_data);
    });
    Base::clear_filtration();
  }

  // Iterator over the simplex tree
  Complex_simplex_iterator get_simplices_iterator_begin() {
    // this specific case works because the range is just a pair of iterators - won't work if range was a vector
    return Base::complex_simplex_range().begin();
  }

  Complex_simplex_iterator get_simplices_iterator_end() {
    // this specific case works because the range is just a pair of iterators - won't work if range was a vector
    return Base::complex_simplex_range().end();
  }

  typename std::vector<Simplex_handle>::const_iterator get_filtration_iterator_begin() {
    // Base::initialize_filtration(); already performed in filtration_simplex_range
    // this specific case works because the range is just a pair of iterators - won't work if range was a vector
    return Base::filtration_simplex_range().begin();
  }

  typename std::vector<Simplex_handle>::const_iterator get_filtration_iterator_end() {
    // this specific case works because the range is just a pair of iterators - won't work if range was a vector
    return Base::filtration_simplex_range().end();
  }

  Skeleton_simplex_iterator get_skeleton_iterator_begin(int dimension) {
    // this specific case works because the range is just a pair of iterators - won't work if range was a vector
    return Base::skeleton_simplex_range(dimension).begin();
  }

  Skeleton_simplex_iterator get_skeleton_iterator_end(int dimension) {
    // this specific case works because the range is just a pair of iterators - won't work if range was a vector
    return Base::skeleton_simplex_range(dimension).end();
  }
};

}  // namespace Gudhi

#endif  // INCLUDE_SIMPLEX_TREE_INTERFACE_H_
