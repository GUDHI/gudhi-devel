/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2014 Inria        """
 *
 *
 *    Modification(s):
 * 		- 2022/11 David Loiseaux / Hannah Schreiber : added multify / flatten to interface standard simplextree.
 *      - YYYY/MM Author: Description of the modification
 */
#ifndef SIMPLEX_TREE_MULTI_H_
#define SIMPLEX_TREE_MULTI_H_

#include <algorithm>
#include <gudhi/Simplex_tree.h>
#include "multi_filtrations/finitely_critical_filtrations.h"
#include "multi_filtrations/line.h"





namespace Gudhi {
/** Model of SimplexTreeOptions.
 *
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */


struct Simplex_tree_options_multidimensional_filtration {
public:
	typedef linear_indexing_tag Indexing_tag;
	typedef int Vertex_handle;
	typedef float value_type;
	using Filtration_value = Gudhi::multi_filtrations::Finitely_critical_multi_filtration<value_type>;
	typedef std::uint32_t Simplex_key;
	static const bool store_key = true;
	static const bool store_filtration = true;
	static const bool contiguous_vertices = false;
	static const bool is_multi_parameter = true;
};



using options_multi = Simplex_tree_options_multidimensional_filtration;
using options_std = Simplex_tree_options_full_featured;
using multi_filtration_type = std::vector<options_multi::value_type>;
using multi_filtration_grid = std::vector<multi_filtration_type>;


template<class options>
Simplex_tree<options>& get_simplextree_from_pointer(const uintptr_t splxptr){ //DANGER
	Simplex_tree<options> &st = *(Gudhi::Simplex_tree<options>*)(splxptr); 
	return st;
}
template<class _options_std, class _options_multi>
void multify(Simplex_tree<_options_std> &st, Simplex_tree<_options_multi> &st_multi, const int dimension){
	if (dimension <= 0)
		{std::cerr << "Empty filtration\n"; throw ;}
	typename _options_multi::Filtration_value f(dimension);
	for (auto &simplex_handle : st.complex_simplex_range()){
		std::vector<int> simplex;
		for (auto vertex : st.simplex_vertex_range(simplex_handle))
			simplex.push_back(vertex);
		f[0] = st.filtration(simplex_handle);
		st_multi.insert_simplex(simplex,f);
	}
}

void multify(const uintptr_t splxptr, const uintptr_t newsplxptr, const int dimension){ //for python
	auto &st = get_simplextree_from_pointer<options_std>(splxptr);
	auto &st_multi = get_simplextree_from_pointer<options_multi>(newsplxptr);
	multify(st, st_multi, dimension);
}

template<class _options_std, class _options_multi>
void flatten(Simplex_tree<_options_std> &st, Simplex_tree<_options_multi> &st_multi, const int dimension = 0){
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		std::vector<int> simplex;
		for (auto vertex : st_multi.simplex_vertex_range(simplex_handle))
			simplex.push_back(vertex);
		typename _options_multi::value_type f = dimension >= 0 ? st_multi.filtration(simplex_handle)[dimension] : 0;
		st.insert_simplex(simplex,f);
	}
}
void flatten(const uintptr_t splxptr, const uintptr_t newsplxptr, const int dimension = 0){ // for python 
	auto &st = get_simplextree_from_pointer<options_std>(newsplxptr);
	auto &st_multi = get_simplextree_from_pointer<options_multi>(splxptr);
	flatten(st, st_multi, dimension);
}

template<class _options_std, class _options_multi>
void flatten_diag(Simplex_tree<_options_std> &st, Simplex_tree<_options_multi> &st_multi, const std::vector<typename _options_multi::value_type> basepoint, int dimension){
	Gudhi::multi_filtrations::Line<typename _options_multi::value_type> l(basepoint);
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		std::vector<int> simplex;
		for (auto vertex : st_multi.simplex_vertex_range(simplex_handle))
			simplex.push_back(vertex);
		
		std::vector<typename _options_multi::value_type> f = st_multi.filtration(simplex_handle);
		if (dimension <0)	 dimension = 0;
		typename _options_multi::value_type new_filtration = l.push_forward(f)[dimension];
		st.insert_simplex(simplex,new_filtration);
	}
}

void flatten_diag(const uintptr_t splxptr, const uintptr_t newsplxptr, const std::vector<options_multi::value_type> basepoint, int dimension){ // for python
	auto &st = get_simplextree_from_pointer<options_std>(newsplxptr);
	auto &st_multi = get_simplextree_from_pointer<options_multi>(splxptr);
	flatten_diag(st,st_multi,basepoint, dimension);
}



template<typename out_type=int>
std::vector<out_type> find_coordinates(const std::vector<options_multi::value_type> &x, const multi_filtration_grid &grid){
	// TODO: optimize with, e.g., dichotomy
	std::vector<out_type> coordinates(grid.size());
	for (int parameter = 0; parameter< (int)grid.size(); parameter++){
		const auto& filtration = grid[parameter];
		const auto& to_project = x[parameter];
		std::vector<options_multi::value_type> distance_vector(filtration.size());
		for (int i = 0; i < (int)filtration.size(); i++){
			distance_vector[i] = std::abs(to_project - filtration[i]);
		}
		coordinates[parameter] = std::distance(distance_vector.begin(), std::min_element(distance_vector.begin(), distance_vector.end()));
	}
	return coordinates;
}

// TODO integer filtrations, does this help with performance ?
// projects filtrations values to the grid. If coordinate_values is set to true, the filtration values are the coordinates of this grid
void squeeze_filtration(uintptr_t splxptr, const multi_filtration_grid &grid, bool coordinate_values=true){
	
	Simplex_tree<options_multi> &st_multi = *(Gudhi::Simplex_tree<options_multi>*)(splxptr);
	int num_parameters = st_multi.get_number_of_parameters();
	if ((int)grid.size() != num_parameters){
		std::cerr << "Bad grid !" << std::endl;
		throw;
	}
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		std::vector<options_multi::value_type> simplex_filtration = st_multi.filtration(simplex_handle);
		if (coordinate_values)
			st_multi.assign_filtration(simplex_handle, find_coordinates<options_multi::value_type>(simplex_filtration, grid));
		else{
			auto coordinates = find_coordinates<int>(simplex_filtration, grid);
			std::vector<options_multi::value_type> squeezed_filtration(num_parameters);
			for (int parameter = 0; parameter < num_parameters; parameter++)
				squeezed_filtration[parameter] = grid[parameter][coordinates[parameter]];
			st_multi.assign_filtration(simplex_handle, squeezed_filtration);
		}
	}
	return;
}
std::vector<multi_filtration_grid> get_filtration_values(const uintptr_t splxptr, const std::vector<int> &degrees){
	Simplex_tree<options_multi> &st_multi = *(Gudhi::Simplex_tree<options_multi>*)(splxptr);
	int num_parameters = st_multi.get_number_of_parameters();
	std::vector<multi_filtration_grid> out(degrees.size(), multi_filtration_grid(num_parameters));
	std::vector<int> degree_index(degrees.size());
	int count = 0;
	for (auto degree : degrees){
		degree_index[degree] = count; count++;
		out[degree_index[degree]].reserve(st_multi.num_simplices());
	}
		
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		const auto filtration = st_multi.filtration(simplex_handle);
		const auto degree = st_multi.dimension(simplex_handle);
		if (std::find(degrees.begin(), degrees.end(), degree) == degrees.end()) continue;
		for (int parameter=0; parameter < num_parameters; parameter++){
			out[degree_index[degree]][parameter].push_back(filtration[parameter]);
		}
	}
	return out;

}

}	// namespace Gudhi

namespace std {

template<>
class numeric_limits<Gudhi::multi_filtration_type>
{
public:
	static Gudhi::multi_filtration_type infinity() throw(){
		return Gudhi::multi_filtration_type(1, std::numeric_limits<Gudhi::Simplex_tree_options_multidimensional_filtration::value_type>::infinity());
	};


	static Gudhi::multi_filtration_type quiet_NaN() throw(){
		return Gudhi::multi_filtration_type(1, numeric_limits<Gudhi::Simplex_tree_options_multidimensional_filtration::value_type>::quiet_NaN());
	};

};

}


#endif  // SIMPLEX_TREE_MULTI_H_
