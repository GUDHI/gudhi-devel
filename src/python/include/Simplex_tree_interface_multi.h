/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *		- 2022/11 David Loiseaux, Hannah Schreiber : adapt for multipersistence. 
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_SIMPLEX_TREE_INTERFACE_H_
#define INCLUDE_SIMPLEX_TREE_INTERFACE_H_

//#include <gudhi/graph_simplicial_complex.h>
//#include <gudhi/distance_functions.h>
//#include <gudhi/Points_off_io.h>
//#include <gudhi/Flag_complex_edge_collapser.h>
#include <gudhi/Simplex_tree.h>
#include "Simplex_tree_multi.h"
// #include <omp.h>
#include "multi_filtrations/finitely_critical_filtrations.h"

#include <iostream>
#include <vector>
#include <utility>  // std::pair
#include <tuple>
#include <iterator>  // for std::distance

namespace Gudhi {

template<typename SimplexTreeOptions = Simplex_tree_options_full_featured>
class Simplex_tree_interface : public Simplex_tree<SimplexTreeOptions> {
 public:
  using Python_filtration_type = std::vector<typename SimplexTreeOptions::value_type>; // TODO : std::conditional
  using Base = Simplex_tree<SimplexTreeOptions>;
  using Filtration_value = typename Base::Filtration_value;
  using Vertex_handle = typename Base::Vertex_handle;
  using Simplex_handle = typename Base::Simplex_handle;
  using Insertion_result = typename std::pair<Simplex_handle, bool>;
  using Simplex = std::vector<Vertex_handle>;
  using Simplex_and_filtration = std::pair<Simplex, Python_filtration_type>;
  using Filtered_simplices = std::vector<Simplex_and_filtration>;
  using Skeleton_simplex_iterator = typename Base::Skeleton_simplex_iterator;
  using Complex_simplex_iterator = typename Base::Complex_simplex_iterator;
  using Extended_filtration_data = typename Base::Extended_filtration_data;
  using Boundary_simplex_iterator = typename Base::Boundary_simplex_iterator;
  typedef bool (*blocker_func_t)(Simplex simplex, void *user_data);
  using euler_chars_type = std::vector<int>;

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

  Python_filtration_type simplex_filtration(const Simplex& simplex) {
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

/*  Simplex_tree_interface* collapse_edges(int nb_collapse_iteration) {*/
/*    using Filtered_edge = std::tuple<Vertex_handle, Vertex_handle, Filtration_value>;*/
/*    std::vector<Filtered_edge> edges;*/
/*    for (Simplex_handle sh : Base::skeleton_simplex_range(1)) {*/
/*      if (Base::dimension(sh) == 1) {*/
/*        typename Base::Simplex_vertex_range rg = Base::simplex_vertex_range(sh);*/
/*        auto vit = rg.begin();*/
/*        Vertex_handle v = *vit;*/
/*        Vertex_handle w = *++vit;*/
/*        edges.emplace_back(v, w, Base::filtration(sh));*/
/*      }*/
/*    }*/

/*    for (int iteration = 0; iteration < nb_collapse_iteration; iteration++) {*/
/*      edges = Gudhi::collapse::flag_complex_collapse_edges(std::move(edges));*/
/*    }*/
/*    Simplex_tree_interface* collapsed_stree_ptr = new Simplex_tree_interface();*/
/*    // Copy the original 0-skeleton*/
/*    for (Simplex_handle sh : Base::skeleton_simplex_range(0)) {*/
/*      collapsed_stree_ptr->insert({*(Base::simplex_vertex_range(sh).begin())}, Base::filtration(sh));*/
/*    }*/
/*    // Insert remaining edges*/
/*    for (auto remaining_edge : edges) {*/
/*      collapsed_stree_ptr->insert({std::get<0>(remaining_edge), std::get<1>(remaining_edge)}, std::get<2>(remaining_edge));*/
/*    }*/
/*    return collapsed_stree_ptr;*/
/*  }*/

  void expansion_with_blockers_callback(int dimension, blocker_func_t user_func, void *user_data) {
	Base::expansion_with_blockers(dimension, [&](Simplex_handle sh){
	  Simplex simplex(Base::simplex_vertex_range(sh).begin(), Base::simplex_vertex_range(sh).end());
	  return user_func(simplex, user_data);
	});
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

  std::pair<Boundary_simplex_iterator, Boundary_simplex_iterator> get_boundary_iterators(const Simplex& simplex) {
	auto bd_sh = Base::find(simplex);
	if (bd_sh == Base::null_simplex())
	  throw std::runtime_error("simplex not found - cannot find boundaries");
	// this specific case works because the range is just a pair of iterators - won't work if range was a vector
	auto boundary_srange = Base::boundary_simplex_range(bd_sh);
	return std::make_pair(boundary_srange.begin(), boundary_srange.end());
  }



// ######################## MULTIPERS STUFF
  void set_keys_to_enumerate(){
	int count = 0;
	for (auto sh : Base::filtration_simplex_range())
		Base::assign_key(sh, count++);
  }

  int get_key(const Simplex& simplex){
	return Base::key(Base::find(simplex));
  }

  void set_key(Simplex& simplex, int key){
	Base::assign_key(Base::find(simplex), key);
	return;
  }
  void fill_lowerstar(std::vector<double> filtration, int axis){
	for (auto &SimplexHandle : Base::complex_simplex_range()){
		std::vector<double> current_birth = Base::filtration(SimplexHandle);
		double to_assign = -1*std::numeric_limits<double>::infinity();
		for (auto vertex : Base::simplex_vertex_range(SimplexHandle)){
			to_assign = std::max(filtration[vertex], to_assign);
		}
		current_birth[axis] = to_assign;
		Base::assign_filtration(SimplexHandle, current_birth);
	}
  }
  using simplices_list = std::vector<std::vector<int>>;
  simplices_list get_simplices_of_dimension(int dimension){
	simplices_list simplex_list;
	simplex_list.reserve(Base::num_simplices());
	for (auto &simplexhandle : Base::skeleton_simplex_range(dimension)){
		if (Base::dimension(simplexhandle) == dimension){
			std::vector<int> simplex;
			simplex.reserve(dimension+1);
			for (int vertex : Base::simplex_vertex_range(simplexhandle))
				simplex.push_back(vertex);
			simplex_list.push_back(simplex);
		}
	}
/*	simplex_list.shrink_to_fit();*/
	return simplex_list;
  }
  using edge_list = std::vector<std::pair<std::pair<int,int>, std::pair<double, double>>>;
  edge_list get_edge_list(){
	edge_list simplex_list;
	simplex_list.reserve(Base::num_simplices());
	for (auto &simplexHandle : Base::skeleton_simplex_range(2)){
		if (Base::dimension(simplexHandle) == 1){
			std::pair<int,int> simplex;
			auto it = Base::simplex_vertex_range(simplexHandle).begin();
			simplex = {*it, *(++it)};
			auto f = Base::filtration(simplexHandle);
			simplex_list.push_back({simplex, {f[0], f[1]}});
		}
	}
/*	simplex_list.shrink_to_fit();*/
	return simplex_list;
  }

  euler_chars_type euler_char(std::vector<std::vector<options_multi::value_type>> &point_list){ // TODO multi-critical 
		const int npts = point_list.size();
		if (npts == 0){
			return {};
		}
		using Gudhi::multi_filtrations::Finitely_critical_multi_filtration;
		
		euler_chars_type out(point_list.size(), 0.);

		// auto is_greater = [nparameters](const point_type &a, const point_type &b){ //french greater
		// 	for (int i = 0; i< nparameters; i++)
		// 		if( a[i] < b[i])
		// 			return false;
		// 	return true;
		// };
// #pragma omp parallel for
		for (int i = 0; i< npts; i++){ // Maybe add a pragma here for parallel
			auto &euler_char_at_point = out[i];
// #pragma omp parallel for reduction(+:euler_char_at_point) // GUDHI : not possible, need a RANDOM ACCESS ITERATOR
			for(const auto &SimplexHandle : Base::complex_simplex_range()){
				// const Finitely_critical_multi_filtration<options_multi::value_type> &pt = *(Finitely_critical_multi_filtration<options_multi::value_type>*)(&point_list[i]);
				options_multi::Filtration_value filtration = Base::filtration(SimplexHandle);
				// const Finitely_critical_multi_filtration<options_multi::value_type> &filtration = *(Finitely_critical_multi_filtration<options_multi::value_type>*)(&filtration_);
				Finitely_critical_multi_filtration<options_multi::value_type> pt(point_list[i]);
				// if (is_greater(pt, filtration)){
				if (filtration <= pt){
					int sign = Base::dimension(SimplexHandle) %2 ? -1 : 1;  
					euler_char_at_point += sign;
				}
			}
		}
		return out;
	}
	void resize_all_filtrations(int num){ //TODO : that is for 1 critical filtrations
		if (num < 0)	return;
		for(const auto &SimplexHandle : Base::complex_simplex_range()){
			std::vector<options_multi::value_type> new_filtration_value = Base::filtration(SimplexHandle);
			new_filtration_value.resize(num);
			Base::assign_filtration(SimplexHandle, new_filtration_value);
		}
	}

	
};


}  // namespace Gudhi

#endif  // INCLUDE_SIMPLEX_TREE_INTERFACE_H_
