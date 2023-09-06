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

#pragma once

//#include <gudhi/graph_simplicial_complex.h>
//#include <gudhi/distance_functions.h>
//#include <gudhi/Points_off_io.h>
//#include <gudhi/Flag_complex_edge_collapser.h>
#include "Simplex_tree_interface.h"
#include "Simplex_tree_multi.h"
// #include <omp.h>
#include "multi_filtrations/finitely_critical_filtrations.h"

#include <iostream>
#include <vector>
#include <utility>  // std::pair
#include <tuple>
#include <iterator>  // for std::distance
#include <span>

namespace Gudhi {

template<typename SimplexTreeOptions = Simplex_tree_options_multidimensional_filtration>
class Simplex_tree_interface_multi : public Simplex_tree_interface<Simplex_tree_options_multidimensional_filtration> {
 public:
  using Python_filtration_type = std::vector<typename SimplexTreeOptions::value_type>; // TODO : std::conditional
  using Base = Simplex_tree<SimplexTreeOptions>;
  using Filtration_value = typename Base::Filtration_value;
  using Vertex_handle = typename Base::Vertex_handle;
  using Simplex_handle = typename Base::Simplex_handle;
  using Insertion_result = typename std::pair<Simplex_handle, bool>;
  using Simplex = std::vector<Vertex_handle>;
  using Simplex_and_filtration = std::pair<Simplex, typename SimplexTreeOptions::value_type*>;
  using Filtered_simplices = std::vector<Simplex_and_filtration>;
  using Skeleton_simplex_iterator = typename Base::Skeleton_simplex_iterator;
  using Complex_simplex_iterator = typename Base::Complex_simplex_iterator;
  using Extended_filtration_data = typename Base::Extended_filtration_data;
  using Boundary_simplex_iterator = typename Base::Boundary_simplex_iterator;
  typedef bool (*blocker_func_t)(Simplex simplex, void *user_data);
  using euler_chars_type = std::vector<int>;

 public:

  Extended_filtration_data efd;
  
//   bool find_simplex(const Simplex& simplex) {
// 	return (Base::find(simplex) != Base::null_simplex());
//   }

  void assign_simplex_filtration(const Simplex& simplex, const Filtration_value& filtration) {
	Base::assign_filtration(Base::find(simplex), filtration);
	Base::clear_filtration();
  }
//   void assign_simplex_filtration(const Simplex& simplex, const Python_filtration_type& filtration) {
// 	Filtration_value& filtration_ = *(Filtration_value*)(&filtration); // Jardinage for no copy. 
// 	Base::assign_filtration(Base::find(simplex), filtration_);
// 	Base::clear_filtration();
//   }

  bool insert(const Simplex& simplex, const Filtration_value& filtration ) {
	Insertion_result result = Base::insert_simplex_and_subfaces(simplex, filtration);
	if (result.first != Base::null_simplex())
	  Base::clear_filtration();
	return (result.second);
  }

  // Do not interface this function, only used in alpha complex interface for complex creation
  bool insert_simplex(const Simplex& simplex, const Filtration_value& filtration ) {
	Insertion_result result = Base::insert_simplex(simplex, filtration);
	return (result.second);
  }
//   bool insert_simplex(const Simplex& simplex, const Python_filtration_type& filtration ) {
// 	Filtration_value& filtration_ = *(Filtration_value*)(&filtration); // Jardinage for no copy. 
// 	Insertion_result result = Base::insert_simplex(simplex, filtration);
// 	return (result.second);
//   }

  // Do not interface this function, only used in interface for complex creation
  bool insert_simplex_and_subfaces(const Simplex& simplex, const Filtration_value& filtration ) {
	Insertion_result result = Base::insert_simplex_and_subfaces(simplex, filtration);
	return (result.second);
  }
//   bool insert_simplex_and_subfaces(const Simplex& simplex, const Python_filtration_type& filtration ) {
// 	Filtration_value& filtration_ = *(Filtration_value*)(&filtration); // Jardinage for no copy. 
// 	Insertion_result result = Base::insert_simplex_and_subfaces(simplex, filtration);
// 	return (result.second);
//   }

  // Do not interface this function, only used in strong witness interface for complex creation
  bool insert_simplex(const std::vector<std::size_t>& simplex, const Filtration_value& filtration ) {
	Insertion_result result = Base::insert_simplex(simplex, filtration);
	return (result.second);
  }

  // Do not interface this function, only used in strong witness interface for complex creation
  bool insert_simplex_and_subfaces(const std::vector<std::size_t>& simplex, const Filtration_value& filtration) {
	Insertion_result result = Base::insert_simplex_and_subfaces(simplex, filtration);
	return (result.second);
  }
   typename SimplexTreeOptions::value_type* simplex_filtration(const Simplex& simplex) {
	auto& filtration = Base::filtration(Base::find(simplex));
	return &filtration[0];
  }


  Simplex_and_filtration get_simplex_and_filtration(Simplex_handle f_simplex) {
	// Simplex simplex;
	// for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
	// //   simplex.insert(simplex.begin(), vertex); // why not push back ?
	// }
	auto it = Base::simplex_vertex_range(f_simplex);
	Simplex simplex(it.begin(), it.end());
	std::reverse(simplex.begin(), simplex.end());
	return std::make_pair(std::move(simplex), &Base::filtration(f_simplex)[0]);
  }

  Filtered_simplices get_star(const Simplex& simplex) {
	Filtered_simplices star;
	for (auto f_simplex : Base::star_simplex_range(Base::find(simplex))) {
	  Simplex simplex_star;
	  for (auto vertex : Base::simplex_vertex_range(f_simplex)) {
		simplex_star.insert(simplex_star.begin(), vertex);
	  }
	  star.push_back(std::make_pair(simplex_star, &Base::filtration(f_simplex)[0]));
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
	  cofaces.push_back(std::make_pair(simplex_coface, &Base::filtration(f_simplex)[0]));
	}
	return cofaces;
  }

  void compute_extended_filtration() {
	throw std::logic_error("Incompatible with multipers");
  }

  Simplex_tree_interface_multi* collapse_edges(int nb_collapse_iteration) {
	throw std::logic_error("Incompatible with multipers");
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

  void set_key(const Simplex& simplex, int key){
	Base::assign_key(Base::find(simplex), key);
	return;
  }
  void fill_lowerstar(const std::vector<options_multi::value_type>& filtration, int axis){
	using value_type=options_multi::value_type;
	for (auto &SimplexHandle : Base::complex_simplex_range()){
		std::vector<value_type> current_birth = Base::filtration(SimplexHandle);
		value_type to_assign = -1*std::numeric_limits<value_type>::infinity();
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
	for (auto simplexhandle : Base::skeleton_simplex_range(dimension)){
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
	for (auto &simplexHandle : Base::skeleton_simplex_range(1)){
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
// DEPRECATED, USE COORDINATE SIMPLEX TREE
//   euler_chars_type euler_char(const std::vector<std::vector<options_multi::value_type>> &point_list){ // TODO multi-critical 
// 		const int npts = point_list.size();
// 		if (npts == 0){
// 			return {};
// 		}
// 		using Gudhi::multi_filtrations::Finitely_critical_multi_filtration;
		
// 		euler_chars_type out(point_list.size(), 0.);

// 		// auto is_greater = [nparameters](const point_type &a, const point_type &b){ //french greater
// 		// 	for (int i = 0; i< nparameters; i++)
// 		// 		if( a[i] < b[i])
// 		// 			return false;
// 		// 	return true;
// 		// };
// // #pragma omp parallel for
// 		for (int i = 0; i< npts; i++){ // Maybe add a pragma here for parallel
// 			auto &euler_char_at_point = out[i];
// // #pragma omp parallel for reduction(+:euler_char_at_point) // GUDHI : not possible, need a RANDOM ACCESS ITERATOR
// 			for(const auto &SimplexHandle : Base::complex_simplex_range()){
// 				// const Finitely_critical_multi_filtration<options_multi::value_type> &pt = *(Finitely_critical_multi_filtration<options_multi::value_type>*)(&point_list[i]);
// 				options_multi::Filtration_value filtration = Base::filtration(SimplexHandle);
// 				// const Finitely_critical_multi_filtration<options_multi::value_type> &filtration = *(Finitely_critical_multi_filtration<options_multi::value_type>*)(&filtration_);
// 				Finitely_critical_multi_filtration<options_multi::value_type> pt(point_list[i]);
// 				// if (is_greater(pt, filtration)){
// 				if (filtration <= pt){
// 					int sign = Base::dimension(SimplexHandle) %2 ? -1 : 1;  
// 					euler_char_at_point += sign;
// 				}
// 			}
// 		}
// 		return out;
// 	}
	void resize_all_filtrations(int num){ //TODO : that is for 1 critical filtrations
		if (num < 0)	return;
		for(const auto &SimplexHandle : Base::complex_simplex_range()){
			std::vector<options_multi::value_type> new_filtration_value = Base::filtration(SimplexHandle);
			new_filtration_value.resize(num);
			Base::assign_filtration(SimplexHandle, new_filtration_value);
		}
	}

	
};


using interface_std = Simplex_tree_interface<Simplex_tree_options_full_featured>;
using interface_multi = Simplex_tree_interface_multi<Simplex_tree_options_multidimensional_filtration>;


void flatten_diag_from_ptr(const uintptr_t splxptr, const uintptr_t newsplxptr, const std::vector<interface_multi::Options::value_type> basepoint, int dimension){ // for python
	auto &st = get_simplextree_from_pointer<interface_std>(newsplxptr);
	auto &st_multi = get_simplextree_from_pointer<interface_multi>(splxptr);
	flatten_diag(st,st_multi,basepoint, dimension);
}
void multify_from_ptr(uintptr_t splxptr,uintptr_t newsplxptr, const int dimension, const multi_filtration_type& default_values){ //for python
	auto &st = get_simplextree_from_pointer<interface_std>(splxptr);
	auto &st_multi = get_simplextree_from_pointer<interface_multi>(newsplxptr);
	multify(st, st_multi, dimension, default_values);
}
void flatten_from_ptr(uintptr_t splxptr, uintptr_t newsplxptr, const int dimension = 0){ // for python 
	auto &st = get_simplextree_from_pointer<interface_std>(newsplxptr);
	auto &st_multi = get_simplextree_from_pointer<interface_multi>(splxptr);
	flatten(st, st_multi, dimension);
}
template<typename ... Args>
void linear_projection_from_ptr(const uintptr_t ptr, const uintptr_t ptr_multi, Args...args){
	auto &st = get_simplextree_from_pointer<interface_std>(ptr);
	auto &st_multi = get_simplextree_from_pointer<interface_multi>(ptr_multi);
	linear_projection(st, st_multi, args...);
}


}  // namespace Gudhi

