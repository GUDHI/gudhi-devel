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





namespace Gudhi::multiparameter {
/** Model of SimplexTreeOptions.
 *
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */


struct Simplex_tree_options_multidimensional_filtration {
public:
	typedef linear_indexing_tag Indexing_tag;
	typedef int Vertex_handle;
	typedef float value_type;
	using Filtration_value = multi_filtrations::Finitely_critical_multi_filtration<value_type>;
	typedef std::uint32_t Simplex_key;
	static const bool store_key = true;
	static const bool store_filtration = true;
	static const bool contiguous_vertices = false;
	static const bool link_nodes_by_label = true;
	static const bool stable_simplex_handles = false;
	static const bool is_multi_parameter = true;
};



using options_multi = Simplex_tree_options_multidimensional_filtration;
using options_std = Simplex_tree_options_full_featured;
using simplextree_std = Simplex_tree<options_std>;
using simplextree_multi = Simplex_tree<options_multi>;

using multi_filtration_type = std::vector<options_multi::value_type>;
using multi_filtration_grid = std::vector<multi_filtration_type>;

// Retrieves a simplextree from a pointer. As the python simplex_tree and simplex_tree_multi don't know each other we have to exchange them using pointers.
template<class simplextreeinterface>
simplextreeinterface& get_simplextree_from_pointer(const uintptr_t splxptr){ //DANGER
	simplextreeinterface &st = *(simplextreeinterface*)(splxptr); 
	return st;
}

// Turns a 1-parameter simplextree into a multiparameter simplextree, and keeps the 1-filtration in the 1st axis 
template<class simplextree_std, class simplextree_multi>
void multify(simplextree_std &st, simplextree_multi &st_multi, const int num_parameters, const typename simplextree_multi::Options::Filtration_value& default_values={}){
	typename simplextree_multi::Options::Filtration_value f(num_parameters);
	for (auto i = 0u; i<std::min(static_cast<unsigned int>(default_values.size()), static_cast<unsigned int>(num_parameters-1));i++)
		f[i+1] = default_values[i];
	std::vector<int> simplex;
	simplex.reserve(st.dimension()+1);
	for (auto &simplex_handle : st.complex_simplex_range()){
		simplex.clear();
		for (auto vertex : st.simplex_vertex_range(simplex_handle))
			simplex.push_back(vertex);

		if (num_parameters > 0)
		f[0] = st.filtration(simplex_handle);
		auto filtration_copy = f;
		st_multi.insert_simplex(simplex,filtration_copy);
		
	}
}



// Turns a multi-parameter simplextree into a 1-parameter simplextree
template<class simplextree_std, class simplextree_multi>
void flatten(simplextree_std &st, simplextree_multi &st_multi, const int dimension = 0){
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		std::vector<int> simplex;
		for (auto vertex : st_multi.simplex_vertex_range(simplex_handle))
			simplex.push_back(vertex);
		typename simplextree_multi::Options::value_type f = dimension >= 0 ? st_multi.filtration(simplex_handle)[dimension] : 0;
		st.insert_simplex(simplex,f);
	}
}

// Applies a linear form (i.e. scalar product with Riesz rpz) to the filtration to flatten a simplextree multi
template<class simplextree_std, class simplextree_multi>
void linear_projection(simplextree_std &st, simplextree_multi &st_multi, const std::vector<typename simplextree_multi::Options::value_type>& linear_form){
	auto sh = st.complex_simplex_range().begin();
	auto sh_multi = st_multi.complex_simplex_range().begin();
	auto end = st.complex_simplex_range().end();
	typename simplextree_multi::Options::Filtration_value multi_filtration;
	for (; sh != end; ++sh, ++sh_multi){
		multi_filtration = st_multi.filtration(*sh_multi);
		auto projected_filtration = multi_filtration.linear_projection(linear_form);
		st.assign_filtration(*sh, projected_filtration);
	}
}

template<class simplextree_std, class simplextree_multi>
void flatten_diag(simplextree_std &st, simplextree_multi &st_multi, const std::vector<typename simplextree_multi::Options::value_type> basepoint, int dimension){
	multi_filtrations::Line<typename simplextree_multi::Options::value_type> l(basepoint);
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		std::vector<int> simplex;
		for (auto vertex : st_multi.simplex_vertex_range(simplex_handle))
			simplex.push_back(vertex);
		
		std::vector<typename simplextree_multi::Options::value_type> f = st_multi.filtration(simplex_handle);
		if (dimension <0)	 dimension = 0;
		typename simplextree_multi::Options::value_type new_filtration = l.push_forward(f)[dimension];
		st.insert_simplex(simplex,new_filtration);
	}
}





/// @brief turns filtration value x into coordinates in the grid
/// @tparam out_type 
/// @param x 
/// @param grid 
/// @return 
template<typename out_type=int, typename vector_like>
inline void find_coordinates(vector_like& x, const multi_filtration_grid &grid){
	// TODO: optimize with, e.g., dichotomy
	
	for (auto parameter = 0u; parameter < grid.size(); parameter++){
		const auto& filtration = grid[parameter]; // assumes its sorted
		const auto to_project = x[parameter];
		if constexpr (std::numeric_limits<typename vector_like::value_type>::has_infinity)
			if (to_project == std::numeric_limits<typename vector_like::value_type>::infinity()){
				x[parameter] = std::numeric_limits<typename vector_like::value_type>::infinity();
				continue;
			}
		if (to_project >= filtration.back()){
			x[parameter] = filtration.size()-1;
			continue;
		} // deals with infinite value at the end of the grid

		unsigned int i = 0;
		while (to_project > filtration[i] && i<filtration.size()) {
			i++;
		}
		if (i==0)
			x[parameter] = 0;
		else if (i < filtration.size()){
			typename vector_like::value_type d1,d2;
			d1 = std::abs(filtration[i-1] - to_project);
			d2 = std::abs(filtration[i] - to_project);
			x[parameter] = d1<d2 ? i-1 : i;
		}
	}
}


// TODO integer filtrations, does this help with performance ?
// projects filtrations values to the grid. If coordinate_values is set to true, the filtration values are the coordinates of this grid
void squeeze_filtration(uintptr_t splxptr, const multi_filtration_grid &grid, bool coordinate_values=true){
	
	Simplex_tree<options_multi> &st_multi = *(Gudhi::Simplex_tree<options_multi>*)(splxptr);
	auto num_parameters = static_cast<unsigned int>(st_multi.get_number_of_parameters());
	if (grid.size() != num_parameters){
		std::cerr << "Bad grid !" << std::endl;
		throw;
	}
	for (const auto &simplex_handle : st_multi.complex_simplex_range()){
		auto& simplex_filtration = st_multi.filtration_mutable(simplex_handle);
		find_coordinates<options_multi::value_type>(simplex_filtration, grid); // turns the simplexfiltration into coords in the grid
		if (!coordinate_values){
			for (auto parameter = 0u; parameter < num_parameters; parameter++)
				simplex_filtration[parameter] = grid[parameter][simplex_filtration[parameter]];
		}
	}
	return;
}

// retrieves the filtration values of a simplextree. Useful to generate a grid.
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
class numeric_limits<Gudhi::multiparameter::multi_filtration_type>
{
public:
	static Gudhi::multiparameter::multi_filtration_type infinity() throw(){
		return Gudhi::multiparameter::multi_filtration_type(1, std::numeric_limits<Gudhi::multiparameter::Simplex_tree_options_multidimensional_filtration::value_type>::infinity());
	};


	static Gudhi::multiparameter::multi_filtration_type quiet_NaN() throw(){
		return Gudhi::multiparameter::multi_filtration_type(1, numeric_limits<Gudhi::multiparameter::Simplex_tree_options_multidimensional_filtration::value_type>::quiet_NaN());
	};

};

}


#endif  // SIMPLEX_TREE_MULTI_H_
