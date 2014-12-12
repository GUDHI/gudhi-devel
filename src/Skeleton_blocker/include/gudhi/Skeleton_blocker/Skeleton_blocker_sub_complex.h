 /*    This file is part of the Gudhi Library. The Gudhi library
  *    (Geometric Understanding in Higher Dimensions) is a generic C++
  *    library for computational topology.
  *
  *    Author(s):       David Salinas
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
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

#ifndef GUDHI_SKELETON_BLOCKER_SUB_COMPLEX_H
#define GUDHI_SKELETON_BLOCKER_SUB_COMPLEX_H

#include "gudhi/Skeleton_blocker_complex.h"
#include "gudhi/Skeleton_blocker/Skeleton_blocker_simplex.h"
#include "gudhi/Utils.h"

namespace Gudhi{

namespace skbl {

/**
 * @brief Simplicial subcomplex of a complex represented by a skeleton/blockers pair.
 *
 * @details Stores a subcomplex of a simplicial complex.
 * To simplify explanations below, we will suppose that :
 * - K is the root simplicial complex
 * - L is a subcomplex of K.
 *
 * One vertex of K may exists in L but with a different address.
 * To be able to locate the vertices in K from vertices of L, the class
 * stores a map 'adresses' between vertices of K and vertices of L.
 *
 * Note that the type for handle of vertices of L is 'Vertex_handle' and
 * the type for handle of vertices of K is 'Root_vertex_handle'.
 *
 * The template ComplexType is type of the root complex. It allows to know
 * if the subcomplex is geometric or not.
 * It has to be either 'Skeleton_blockers_complex' or 'Skeleton_blockers_geometric_complex'.
 *
 */
template<typename ComplexType>
class Skeleton_blocker_sub_complex : public ComplexType
{

protected:

	template<class T> friend class Skeleton_blocker_link_complex;

	typedef typename ComplexType::Graph Graph;
	typedef typename ComplexType::Edge_handle Edge_handle;

	typedef typename ComplexType::boost_vertex_handle boost_vertex_handle;


public:
	using ComplexType::add_vertex;
	using ComplexType::add_edge;
	using ComplexType::add_blocker;

	typedef typename ComplexType::Vertex_handle Vertex_handle;
	typedef typename ComplexType::Root_vertex_handle Root_vertex_handle;
	typedef typename ComplexType::Simplex_handle Simplex_handle;
	typedef typename ComplexType::Root_simplex_handle Root_simplex_handle;


protected:


	///**
	//* @brief Returns true iff the simplex formed by all vertices contained in 'addresses_sigma_in_link'
	//* but 'vertex_to_be_ignored' is in 'link'
	//*/
	/*
			template<typename T> friend bool
				proper_face_in_union(
				Skeleton_blocker_sub_complex<T> & link,
				std::vector<boost::optional<typename T::Vertex_handle> > & addresses_sigma_in_link,
				int vertex_to_be_ignored);*/

	/**
	 * @brief Determines whether all proper faces of simplex 'sigma' belong to 'link1' \cup 'link2'
	 * where 'link1' and 'link2' are subcomplexes of the same complex of type ComplexType
	 */
	//	template<typename T> friend	bool
	//	proper_faces_in_union(Skeleton_blocker_simplex<typename T::Root_vertex_handle> & sigma, Skeleton_blocker_sub_complex<T> & link1, Skeleton_blocker_sub_complex<T> & link2){

	//template<typename T> friend bool
	//proper_faces_in_union(Skeleton_blocker_simplex<typename T::Root_vertex_handle> & sigma, Skeleton_blocker_sub_complex<T> & link1, Skeleton_blocker_sub_complex<T> & link2);

	typedef std::map<Root_vertex_handle, Vertex_handle> IdAddressMap;
	typedef typename IdAddressMap::value_type AddressPair;
	typedef typename IdAddressMap::iterator IdAddressMapIterator;
	typedef typename IdAddressMap::const_iterator IdAddressMapConstIterator;
	std::map<Root_vertex_handle, Vertex_handle> adresses;


public:
	/**
	 * Add a vertex 'global' of K to L. When added to L, this vertex will receive
	 * another number, addresses(global), its local adress.
	 * return the adress where the vertex lay on L.
	 * The vertex corresponding to 'global' must not be already present
	 * in the complex.
	 */
	Vertex_handle add_vertex(Root_vertex_handle global){
		assert(!this->contains_vertex(global));
		Vertex_handle address(boost::add_vertex(this->skeleton));
		this->num_vertices_++;
		(*this)[address].activate();
		(*this)[address].set_id(global);
		adresses.insert(AddressPair(global, address));
		this->degree_.push_back(0);
		return address;
	}


	/**
	 * Add an edge (v1_root,v2_root) to the sub-complex.
	 * It assumes that both vertices corresponding to v1_root and v2_root are present
	 * in the sub-complex.
	 */
	void add_edge(Root_vertex_handle v1_root, Root_vertex_handle v2_root){
		auto v1_sub(this->get_address(v1_root));
		auto v2_sub(this->get_address(v2_root));
		assert(v1_sub && v2_sub);
		this->ComplexType::add_edge(*v1_sub, *v2_sub);
	}

	/**
	 * Add a blocker to the sub-complex.
	 * It assumes that all vertices of blocker_root are present
	 * in the sub-complex.
	 */
	void add_blocker(const Root_simplex_handle& blocker_root){
		auto blocker_sub = this->get_address(blocker_root);
		assert(blocker_sub);
		this->add_blocker(new Simplex_handle(*blocker_sub));
	}



public:

	/**
	 * Constructs the restricted complex of 'parent_complex' to
	 * vertices of 'simplex'.
	 */
	friend void make_restricted_complex(const ComplexType & parent_complex, const Simplex_handle& simplex, Skeleton_blocker_sub_complex & result){
		result.clear();
		// add vertices to the sub complex
		for (auto x : simplex){
			assert(parent_complex.contains_vertex(x));
			auto x_local = result.add_vertex(parent_complex[x].get_id());
			result[x_local] = parent_complex[x];
		}

		// add edges to the sub complex
		for (auto x : simplex){
			// x_neigh is the neighbor of x intersected with vertices_simplex
			Simplex_handle x_neigh;
			parent_complex.add_neighbours(x, x_neigh, true);
			x_neigh.intersection(simplex);
			for (auto y : x_neigh){
				result.add_edge(parent_complex[x].get_id(), parent_complex[y].get_id());
			}
		}

		// add blockers to the sub complex
		for (auto blocker : parent_complex.const_blocker_range()){
			// check if it is the first time we encounter the blocker
			if (simplex.contains(*blocker)){
				Root_simplex_handle blocker_root(parent_complex.get_id(*(blocker)));
				Simplex_handle blocker_restr(*result.get_simplex_address(blocker_root));
				result.add_blocker(new Simplex_handle(blocker_restr));
			}
		}
	}



	void clear(){
		adresses.clear();
		ComplexType::clear();
	}

	/**
	 * Compute the local vertex in L corresponding to the vertex global in K.
	 * runs in O(log n) if n = num_vertices()
	 */
	boost::optional<Vertex_handle> get_address(Root_vertex_handle global) const{
		boost::optional<Vertex_handle> res;
		IdAddressMapConstIterator it = adresses.find(global);
		if (it == adresses.end()) res.reset();
		else  res = (*it).second;
		return res;
	}

	//	/**
	//	 * Allocates a simplex in L corresponding to the simplex s in K
	//	 * with its local adresses and returns an AddressSimplex.
	//	 */
	//	boost::optional<Simplex_handle> get_address(const Root_simplex_handle & s) const;


//private:
	/**
	 *  same as get_address except that it will return a simplex in any case.
	 *  The vertices that were not found are not added.
	 */
	// @remark should be private but problem with VS
	std::vector<boost::optional<Vertex_handle> > get_addresses(const Root_simplex_handle & s) const{
		std::vector<boost::optional<Vertex_handle> > res;
		for (auto i : s)
		{
			res.push_back(get_address(i));
		}
		return res;
	}

};


/**
 * @remark remarque perte de temps a creer un nouveau simplexe a chaque fois
 * alors qu'on pourrait utiliser a la place de 'addresses_sigma_in_link'
 * un simplex avec des valeurs spéciales ComplexDS::null_vertex par exemple
 * pour indiquer qu'un vertex n'appartient pas au complex
 */
template<typename ComplexType>
bool proper_face_in_union(
		Skeleton_blocker_sub_complex<ComplexType> & link,
		std::vector<boost::optional<typename ComplexType::Vertex_handle> > & addresses_sigma_in_link,
		int vertex_to_be_ignored)
{
	// we test that all vertices of 'addresses_sigma_in_link' but 'vertex_to_be_ignored'
	// are in link1 if it is the case we construct the corresponding simplex
	bool vertices_sigma_are_in_link = true;
	typename ComplexType::Simplex_handle sigma_in_link;
	for (int i = 0; i < addresses_sigma_in_link.size(); ++i){
		if (i != vertex_to_be_ignored){
			if (!addresses_sigma_in_link[i]){
				vertices_sigma_are_in_link = false;
				break;
			}
			else sigma_in_link.add_vertex(*addresses_sigma_in_link[i]);
		}
	}
	// If one of vertices of the simplex is not in the complex then it returns false
	// Otherwise, it tests if the simplex is in the complex
	return vertices_sigma_are_in_link && link.contains(sigma_in_link);
}

/*
		template<typename ComplexType>
		bool
			proper_faces_in_union(Skeleton_blocker_simplex<typename ComplexType::Root_vertex_handle> & sigma, Skeleton_blocker_sub_complex<ComplexType> & link1, Skeleton_blocker_sub_complex<ComplexType> & link2)
		{
				typedef typename ComplexType::Vertex_handle  Vertex_handle;
				std::vector<boost::optional<Vertex_handle> > addresses_sigma_in_link1 = link1.get_addresses(sigma);
				std::vector<boost::optional<Vertex_handle> > addresses_sigma_in_link2 = link2.get_addresses(sigma);

				for (int current_index = 0; current_index < addresses_sigma_in_link1.size(); ++current_index)
				{

					if (!proper_face_in_union(link1, addresses_sigma_in_link1, current_index)
						&& !proper_face_in_union(link2, addresses_sigma_in_link2, current_index)){
						return false;
					}
				}
				return true;
			}*/



// Remark: this function should be friend in order to leave get_adresses private
// however doing so seemes currently not possible due to a visual studio bug c2668 
// "the compiler does not support partial ordering of template functions as specified in the C++ Standard"
// http://www.serkey.com/error-c2668-ambiguous-call-to-overloaded-function-bb45ft.html
template<typename ComplexType>
bool
proper_faces_in_union(Skeleton_blocker_simplex<typename ComplexType::Root_vertex_handle> & sigma, Skeleton_blocker_sub_complex<ComplexType> & link1, Skeleton_blocker_sub_complex<ComplexType> & link2)
{
	typedef typename ComplexType::Vertex_handle  Vertex_handle;
	std::vector<boost::optional<Vertex_handle> > addresses_sigma_in_link1 = link1.get_addresses(sigma);
	std::vector<boost::optional<Vertex_handle> > addresses_sigma_in_link2 = link2.get_addresses(sigma);

	for (int current_index = 0; current_index < addresses_sigma_in_link1.size(); ++current_index)
	{

		if (!proper_face_in_union(link1, addresses_sigma_in_link1, current_index)
				&& !proper_face_in_union(link2, addresses_sigma_in_link2, current_index)){
			return false;
		}
	}
	return true;
}

} // namespace skbl

}  // namespace GUDHI


#endif

