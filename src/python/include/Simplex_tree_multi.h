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
#include "utils/box.h"
#include "utils/line.h"



namespace Gudhi {
/** Model of SimplexTreeOptions.
 *
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */

using value_type = double;
struct Simplex_tree_options_multidimensional_filtration {
public:
	typedef linear_indexing_tag Indexing_tag;
	typedef int Vertex_handle;
	typedef std::vector<value_type> Filtration_value;
	typedef std::uint32_t Simplex_key;
	static const bool store_key = true;
	static const bool store_filtration = true;
	static const bool contiguous_vertices = false;
};

using option_multi = Simplex_tree_options_multidimensional_filtration;
using option_std = Simplex_tree_options_full_featured;
bool operator<(const std::vector<value_type>& v1, const std::vector<value_type>& v2)
{
	bool isSame = true;
	if (v1.size() != v2.size()) isSame = false;
	for (unsigned int i = 0; i < std::min(v1.size(), v2.size()); ++i){
		if (v1[i] > v2[i]) return false;
		if (isSame && v1[i] != v2[i]) isSame = false;
	}
	if (isSame) return false;
	return true;
}

void multify(const uintptr_t splxptr, const uintptr_t newsplxptr, const unsigned int dimension){
	Simplex_tree<option_std> &st = *(Gudhi::Simplex_tree<option_std>*)(splxptr);
	Simplex_tree<option_multi> &st_multi = *(Gudhi::Simplex_tree<option_multi>*)(newsplxptr);;
	if (dimension <= 0)
		{std::cout << "Empty filtration\n"; return ;}
	std::vector<value_type> f(dimension);
	for (auto &simplex_handle : st.complex_simplex_range()){
		std::vector<int> simplex;
		for (auto vertex : st.simplex_vertex_range(simplex_handle))
			simplex.push_back(vertex);
		f[0] = st.filtration(simplex_handle);
		st_multi.insert_simplex(simplex,f);
	}
}
void flatten(const uintptr_t splxptr, const uintptr_t newsplxptr, const unsigned int dimension = 0){
	Simplex_tree<option_std> &st = *(Gudhi::Simplex_tree<option_std>*)(newsplxptr);
	Simplex_tree<option_multi> &st_multi = *(Gudhi::Simplex_tree<option_multi>*)(splxptr);

	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		std::vector<int> simplex;
		for (auto vertex : st_multi.simplex_vertex_range(simplex_handle))
			simplex.push_back(vertex);
		value_type f = st_multi.filtration(simplex_handle)[dimension];
		st.insert_simplex(simplex,f);
	}
}

void flatten_diag(const uintptr_t splxptr, const uintptr_t newsplxptr, const std::vector<value_type> basepoint, int dimension){
	Simplex_tree<option_std> &st = *(Gudhi::Simplex_tree<option_std>*)(newsplxptr);
	Simplex_tree<option_multi> &st_multi = *(Gudhi::Simplex_tree<option_multi>*)(splxptr);
	utils::Line l(basepoint);
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		std::vector<int> simplex;
		for (auto vertex : st_multi.simplex_vertex_range(simplex_handle))
			simplex.push_back(vertex);
		
		std::vector<value_type> f = st_multi.filtration(simplex_handle);
		if (dimension <0)	 dimension = 0;
		value_type new_filtration = l.push_forward(f)[dimension];
		st.insert_simplex(simplex,new_filtration);
	}


}


using filtration_grid = std::vector<std::vector<value_type>>;
using grid_value = std::vector<value_type>;
template<typename out_type=int>
std::vector<out_type> find_coordinates(const grid_value &x, const filtration_grid &grid){
	// TODO: optimize with dichotomy
	std::vector<out_type> coordinates(grid.size());
	for (unsigned int parameter = 0; parameter< grid.size(); parameter++){
		const auto& filtration = grid[parameter];
		const auto& to_project = x[parameter];
		grid_value distance_vector(filtration.size());
		for (unsigned int i = 0; i < filtration.size(); i++){
			distance_vector[i] = std::abs(to_project - filtration[i]);
		}
		coordinates[parameter] = std::distance(distance_vector.begin(), std::min_element(distance_vector.begin(), distance_vector.end()));
	}
	return coordinates;
}

// TODO integer filtrations, does this help with performance ?
// projects filtrations values to the grid. If coordinate_values is set to true, the filtration values are the coordinates of this grid
void squeeze_filtration(uintptr_t splxptr, const std::vector<std::vector<value_type>> &grid, bool coordinate_values=true){
	
	Simplex_tree<option_multi> &st_multi = *(Gudhi::Simplex_tree<option_multi>*)(splxptr);
	unsigned int num_parameters = st_multi.get_number_of_parameters();
	if (grid.size() != num_parameters){
		std::cerr << "Bad grid !" << std::endl;
		return;
	}
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		std::vector<value_type> simplex_filtration = st_multi.filtration(simplex_handle);
		if (coordinate_values)
			st_multi.assign_filtration(simplex_handle, find_coordinates<value_type>(simplex_filtration, grid));
		else{
			auto coordinates = find_coordinates<int>(simplex_filtration, grid);
			grid_value squeezed_filtration(num_parameters);
			for (unsigned int parameter = 0; parameter < num_parameters; parameter++)
				squeezed_filtration[parameter] = grid[parameter][coordinates[parameter]];
			st_multi.assign_filtration(simplex_handle, squeezed_filtration);
		}
	}
	return;
}
std::vector<std::vector<value_type>> get_filtration_values(const uintptr_t splxptr){
	Simplex_tree<option_multi> &st_multi = *(Gudhi::Simplex_tree<option_multi>*)(splxptr);
	unsigned int num_parameters = st_multi.get_number_of_parameters();
	std::vector<std::vector<value_type>> out(num_parameters, std::vector<value_type>(st_multi.num_simplices()));
	unsigned int count = 0;
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		const auto filtration = st_multi.filtration(simplex_handle);
		for (unsigned int parameter=0; parameter < num_parameters; parameter++){
			out[parameter][count] = filtration[parameter];
		}
		count++;
	}
	return out;

}

}	// namespace Gudhi

namespace std {

template<>
class numeric_limits<std::vector<Gudhi::value_type> >
{
public:
	static std::vector<Gudhi::value_type> infinity() throw(){
		return std::vector<Gudhi::value_type>(1, numeric_limits<Gudhi::value_type>::infinity());
	};


	static std::vector<Gudhi::value_type> quiet_NaN() throw(){
		return std::vector<Gudhi::value_type>(1, numeric_limits<Gudhi::value_type>::quiet_NaN());
	};

};

}	// namespace std





#endif  // SIMPLEX_TREE_MULTI_H_
