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
#ifndef GUDHI_SKELETON_BLOCKERS_SIMPLIFIABLE_COMPLEX_H_
#define GUDHI_SKELETON_BLOCKERS_SIMPLIFIABLE_COMPLEX_H_

#include "gudhi/Skeleton_blocker/Skeleton_blocker_sub_complex.h"

namespace Gudhi{

namespace skbl {

/**
 *  \brief Class that allows simplification operation on a simplicial complex represented
 *  by a skeleton/blockers pair.
 */
template<typename SkeletonBlockerDS>
class Skeleton_blocker_simplifiable_complex : public Skeleton_blocker_complex<SkeletonBlockerDS>
{
	template<class ComplexType> friend class Skeleton_blocker_sub_complex;

public:

	typedef Skeleton_blocker_complex<SkeletonBlockerDS> SkeletonBlockerComplex;

	typedef typename SkeletonBlockerComplex::Graph_edge Graph_edge;

	typedef typename SkeletonBlockerComplex::boost_adjacency_iterator boost_adjacency_iterator;
	typedef typename SkeletonBlockerComplex::Edge_handle Edge_handle;
	typedef typename SkeletonBlockerComplex::boost_vertex_handle boost_vertex_handle;
	typedef typename SkeletonBlockerComplex::Vertex_handle Vertex_handle;
	typedef typename SkeletonBlockerComplex::Root_vertex_handle Root_vertex_handle;
	typedef typename SkeletonBlockerComplex::Simplex_handle Simplex_handle;
	typedef typename SkeletonBlockerComplex::Root_simplex_handle Root_simplex_handle;
	typedef typename SkeletonBlockerComplex::Blocker_handle Blocker_handle;


	typedef typename SkeletonBlockerComplex::Root_simplex_iterator Root_simplex_iterator;
	typedef typename SkeletonBlockerComplex::Simplex_handle_iterator Simplex_handle_iterator;
	typedef typename SkeletonBlockerComplex::BlockerMap BlockerMap;
	typedef typename SkeletonBlockerComplex::BlockerPair BlockerPair;
	typedef typename SkeletonBlockerComplex::BlockerMapIterator BlockerMapIterator;
	typedef typename SkeletonBlockerComplex::BlockerMapConstIterator BlockerMapConstIterator;

	typedef typename SkeletonBlockerComplex::Visitor Visitor;


	/** @name Constructors / Destructors / Initialization
	 */
	//@{
	Skeleton_blocker_simplifiable_complex(int num_vertices_ = 0,Visitor* visitor_=NULL):
		Skeleton_blocker_complex<SkeletonBlockerDS>(num_vertices_,visitor_){	}


	/**
	 * @brief Constructor with a list of simplices
	 * @details The list of simplices must be the list
	 * of simplices of a simplicial complex, sorted with increasing dimension.
	 * todo take iterator instead
	 */
	Skeleton_blocker_simplifiable_complex(std::list<Simplex_handle>& simplices,Visitor* visitor_=NULL):
		Skeleton_blocker_complex<SkeletonBlockerDS>(simplices,visitor_)
		{}



	virtual ~Skeleton_blocker_simplifiable_complex(){
	}

	//@}

	/**
	 * Returns true iff the blocker 'sigma' is popable.
	 * To define popable, let us call 'L' the complex that
	 * consists in the current complex without the blocker 'sigma'.
	 * A blocker 'sigma' is then "popable" if the link of 'sigma'
	 * in L is reducible.
	 *
	 */
	virtual bool is_popable_blocker(Blocker_handle sigma) const{
		assert(this->contains_blocker(*sigma));
		Skeleton_blocker_link_complex<Skeleton_blocker_simplifiable_complex> link_blocker_sigma;
		build_link_of_blocker(*this,*sigma,link_blocker_sigma);

		bool res = link_blocker_sigma.is_contractible()==CONTRACTIBLE;
		return res;
	}




private:
	/**
	 * @returns the list of blockers of the simplex
	 *
	 * @todo a enlever et faire un iterateur sur tous les blockers a la place
	 */
	std::list<Blocker_handle> get_blockers(){
		std::list<Blocker_handle> res;
		for (auto blocker : this->blocker_range()){
			res.push_back(blocker);
		}
		return res;
	}



public:

	/**
	 * Removes all the popable blockers of the complex and delete them.
	 * @returns the number of popable blockers deleted
	 */
	void remove_popable_blockers(){
		std::list<Vertex_handle> vertex_to_check;
		for(auto v : this->vertex_range())
			vertex_to_check.push_front(v);

		while(!vertex_to_check.empty()){
			Vertex_handle v = vertex_to_check.front();
			vertex_to_check.pop_front();

			bool blocker_popable_found=true;
			while (blocker_popable_found){
				blocker_popable_found = false;
				for(auto block : this->blocker_range(v)){
					if (this->is_popable_blocker(block)) {
						for(Vertex_handle nv : *block)
							if(nv!=v) vertex_to_check.push_back(nv);
						this->delete_blocker(block);
						blocker_popable_found = true;
						break;
					}
				}
			}
		}
	}


	/**
	 * Removes all the popable blockers of the complex passing through v and delete them.
	 */
	void remove_popable_blockers(Vertex_handle v){
		bool blocker_popable_found=true;
		while (blocker_popable_found){
			blocker_popable_found = false;
			for(auto block : this->blocker_range(v)){
				if (is_popable_blocker(block)) {
					this->delete_blocker(block);
					blocker_popable_found = true;
				}
			}
		}
	}


	/**
	 * Remove the star of the vertex 'v'
	 */
	void remove_star(Vertex_handle v){
		// we remove the blockers that are not consistent anymore

		update_blockers_after_remove_star_of_vertex_or_edge(v);

		while (this->degree(v) > 0)
		{
			Vertex_handle w( * (adjacent_vertices(v.vertex, this->skeleton).first));
			this->remove_edge(v,w);
		}
		this->remove_vertex(v);
	}

private:
	/**
	 * after removing the star of a simplex, blockers sigma that contains this simplex must be removed.
	 * Furthermore, all simplices tau of the form sigma \setminus simplex_to_be_removed must be added
	 * whenever the dimension of tau is at least 2.
	 */
	void update_blockers_after_remove_star_of_vertex_or_edge(const Simplex_handle& simplex_to_be_removed){
		std::list <Blocker_handle> blockers_to_update;
		if(simplex_to_be_removed.empty()) return;

		auto v0 = simplex_to_be_removed.first_vertex();
		for (auto blocker : this->blocker_range(v0)){
			if(blocker->contains(simplex_to_be_removed))
				blockers_to_update.push_back(blocker);
		}

		for(auto blocker_to_update : blockers_to_update){
			Simplex_handle sub_blocker_to_be_added;
			bool sub_blocker_need_to_be_added =
					(blocker_to_update->dimension()-simplex_to_be_removed.dimension()) >= 2;
			if(sub_blocker_need_to_be_added){
				sub_blocker_to_be_added = *blocker_to_update;
				sub_blocker_to_be_added.difference(simplex_to_be_removed);
			}
			this->delete_blocker(blocker_to_update);
			if(sub_blocker_need_to_be_added)
				this->add_blocker(sub_blocker_to_be_added);
		}
	}



public:
	/**
	 * Remove the star of the edge connecting vertices a and b.
	 * @returns the number of blocker that have been removed
	 */
	void remove_star(Vertex_handle a, Vertex_handle b){
		update_blockers_after_remove_star_of_vertex_or_edge(Simplex_handle(a,b));
		// we remove the edge
		this->remove_edge(a,b);
	}


	/**
	 * Remove the star of the edge 'e'.
	 */
	void remove_star(Edge_handle e){
		return remove_star(this->first_vertex(e),this->second_vertex(e));
	}

	/**
	 * Remove the star of the simplex 'sigma' which needs to belong to the complex
	 */
	void remove_star(const Simplex_handle& sigma){
		assert(this->contains(sigma));
		if (sigma.dimension()==0)
			remove_star(sigma.first_vertex());
		else
			if (sigma.dimension()==1)
				remove_star(sigma.first_vertex(),sigma.last_vertex());
			else
				update_blockers_after_remove_star_of_simplex(sigma);
	}

private:
	void update_blockers_after_remove_star_of_simplex(const Simplex_handle& sigma){
		std::list <Blocker_handle> blockers_to_remove;
		for (auto blocker : this->blocker_range(sigma.first_vertex())){
			if(blocker->contains(sigma))
				blockers_to_remove.push_back(blocker);
		}
		for(auto blocker_to_update : blockers_to_remove)
			this->delete_blocker(blocker_to_update);
		this->add_blocker(sigma);
	}

public:

	enum simplifiable_status{ NOT_HOMOTOPY_EQ,MAYBE_HOMOTOPY_EQ,HOMOTOPY_EQ};
	simplifiable_status is_remove_star_homotopy_preserving(const Simplex_handle& simplex){
		return MAYBE_HOMOTOPY_EQ;
	}



	enum contractible_status{ NOT_CONTRACTIBLE,MAYBE_CONTRACTIBLE,CONTRACTIBLE};
	/**
	 * @brief %Test if the complex is reducible using a strategy defined in the class
	 * (by default it tests if the complex is a cone)
	 * @details Note that NO could be returned if some invariant ensures that the complex
	 * is not a point (for instance if the euler characteristic is different from 1).
	 * This function will surely have to return MAYBE in some case because the
	 * associated problem is undecidable but it in practice, it can often
	 * be solved with the help of geometry.
	 */
	virtual contractible_status is_contractible() const{
		if (this->is_cone()) return CONTRACTIBLE;
		else return MAYBE_CONTRACTIBLE;
		//		return this->is_cone();
	}


	/** @Edge contraction operations
	 */
	//@{



	/**
	 * @return If ignore_popable_blockers is true
	 * then the result is true iff the link condition at edge ab is satisfied
	 * or equivalently iff no blocker contains ab.
	 * If ignore_popable_blockers is false then the
	 * result is true iff all blocker containing ab are popable.
	 */
	bool link_condition(Vertex_handle a, Vertex_handle b,bool ignore_popable_blockers = false) const{
		for (auto blocker : this->const_blocker_range(a))
			if ( blocker->contains(b) ){
				// false if ignore_popable_blockers is false
				// otherwise the blocker has to be popable
				return ignore_popable_blockers && is_popable_blocker(blocker);
			}
		return true;
	}

	/**
	 * @return If ignore_popable_blockers is true
	 * then the result is true iff the link condition at edge ab is satisfied
	 * or equivalently iff no blocker contains ab.
	 * If ignore_popable_blockers is false then the
	 * result is true iff all blocker containing ab are popable.
	 */
	bool link_condition(Edge_handle & e,bool ignore_popable_blockers = false) const{
		const Graph_edge& edge = (*this)[e];
		assert(this->get_address(edge.first()));
		assert(this->get_address(edge.second()));
		Vertex_handle a(*this->get_address(edge.first()));
		Vertex_handle b(*this->get_address(edge.second()));
		return link_condition(a,b,ignore_popable_blockers);
	}


protected:
	/**
	 * Compute simplices beta such that a.beta is an order 0 blocker
	 * that may be used to construct a new blocker after contracting ab.
	 * Suppose that the link condition Link(ab) = Link(a) inter Link(b)
	 * is satisfied.
	 */
	void tip_blockers(Vertex_handle a, Vertex_handle b, std::vector<Simplex_handle> & buffer) const{
		for (auto const & blocker : this->const_blocker_range(a))
		{
			Simplex_handle beta = (*blocker);
			beta.remove_vertex(a);
			buffer.push_back(beta);
		}

		Simplex_handle n;
		this->add_neighbours(b,n);
		this->remove_neighbours(a,n);
		n.remove_vertex(a);


		for (Vertex_handle y : n)
		{
			Simplex_handle beta;
			beta.add_vertex( y );
			buffer.push_back(beta);
		}
	}


private:

	/**
	 * @brief "Replace" the edge 'bx' by the edge 'ax'.
	 * Assume that the edge 'bx' was present whereas 'ax' was not.
	 * Precisely, it does not replace edges, but remove 'bx' and then add 'ax'.
	 * The visitor 'on_swaped_edge' is called just after edge 'ax' had been added
	 * and just before edge 'bx' had been removed. That way, it can
	 * eventually access to information of 'ax'.
	 */
	void swap_edge(Vertex_handle a,Vertex_handle b,Vertex_handle x){
		this->add_edge(a,x);
		if (this->visitor) this->visitor->on_swaped_edge(a,b,x);
		this->remove_edge(b,x);
	}


private:
	/**
	 * @brief removes all blockers passing through the edge 'ab'
	 */
	void delete_blockers_around_vertex(Vertex_handle v){
		std::list <Blocker_handle> blockers_to_delete;
		for (auto blocker : this->blocker_range(v)){
			blockers_to_delete.push_back(blocker);
		}
		while (!blockers_to_delete.empty()){
			this->remove_blocker(blockers_to_delete.back());
			blockers_to_delete.pop_back();
		}

	}
	/**
	 * @brief removes all blockers passing through the edge 'ab'
	 */
	void delete_blockers_around_edge(Vertex_handle a, Vertex_handle b){
		std::list<Blocker_handle> blocker_to_delete;
		for (auto blocker : this->blocker_range(a))
			if (blocker->contains(b)) blocker_to_delete.push_back(blocker);
		while (!blocker_to_delete.empty())
		{
			this->delete_blocker(blocker_to_delete.back());
			blocker_to_delete.pop_back();
		}
	}

public:

	/**
	 * Contracts the edge.
	 * @remark If the link condition Link(ab) = Link(a) inter Link(b) is not satisfied,
	 * it removes first all blockers passing through 'ab'
	 */
	void contract_edge(Edge_handle edge){
		contract_edge(this->first_vertex(edge),this->second_vertex(edge));
	}


	/**
	 * Contracts the edge connecting vertices a and b.
	 * @remark If the link condition Link(ab) = Link(a) inter Link(b) is not satisfied,
	 * it removes first all blockers passing through 'ab'
	 */
	void contract_edge(Vertex_handle a, Vertex_handle b){
		assert(this->contains_vertex(a));
		assert(this->contains_vertex(b));
		assert(this->contains_edge(a,b));

		// if some blockers passes through 'ab', we remove them.
		if (!link_condition(a,b))
			delete_blockers_around_edge(a,b);

		std::set<Simplex_handle> blockers_to_add;

		get_blockers_to_be_added_after_contraction(a,b,blockers_to_add);

		delete_blockers_around_vertices(a,b);

		update_edges_after_contraction(a,b);

		this->remove_vertex(b);

		notify_changed_edges(a);

		for(auto block : blockers_to_add)
			this->add_blocker(block);

		assert(this->contains_vertex(a));
		assert(!this->contains_vertex(b));
	}


private:

	void get_blockers_to_be_added_after_contraction(Vertex_handle a,Vertex_handle b,std::set<Simplex_handle>& blockers_to_add){
		blockers_to_add.clear();

		typedef Skeleton_blocker_link_complex<Skeleton_blocker_complex<SkeletonBlockerDS> > LinkComplexType;

		LinkComplexType link_a(*this,a);
		LinkComplexType link_b(*this,b);

		std::vector<Simplex_handle> vector_alpha, vector_beta;

		tip_blockers(a,b,vector_alpha);
		tip_blockers(b,a,vector_beta);

		for (auto alpha = vector_alpha.begin(); alpha != vector_alpha.end(); ++alpha){
			for (auto beta = vector_beta.begin(); beta != vector_beta.end(); ++beta)
			{
				Simplex_handle sigma = *alpha; sigma.union_vertices(*beta);
				Root_simplex_handle sigma_id = this->get_id(sigma);
				if ( this->contains(sigma) &&
						proper_faces_in_union<SkeletonBlockerComplex>(sigma_id,link_a,link_b))
				{
					//					Blocker_handle blocker = new Simplex_handle(sigma);
					sigma.add_vertex(a);
					blockers_to_add.insert(sigma);
				}
			}
		}
	}


	/**
	 * delete all blockers that passes through a or b
	 */
	void delete_blockers_around_vertices(Vertex_handle a,Vertex_handle b){
		std::vector<Blocker_handle> blocker_to_delete;
		for(auto bl : this->blocker_range(a))
			blocker_to_delete.push_back(bl);
		for(auto bl : this->blocker_range(b))
			blocker_to_delete.push_back(bl);
		while (!blocker_to_delete.empty())
		{
			this->delete_blocker(blocker_to_delete.back());
			blocker_to_delete.pop_back();
		}
	}


	void update_edges_after_contraction(Vertex_handle a,Vertex_handle b){
		// We update the set of edges
		this->remove_edge(a,b);

		// For all edges {b,x} incident to b,
		// we remove {b,x} and add {a,x} if not already there.
		while (this->degree(b)> 0)
		{
			Vertex_handle x(*(adjacent_vertices(b.vertex, this->skeleton).first));
			if(!this->contains_edge(a,x))
				// we 'replace' the edge 'bx' by the edge 'ax'
				this->swap_edge(a,b,x);
			else
				this->remove_edge(b,x);
		}
	}

	void notify_changed_edges(Vertex_handle a){
		// We notify the visitor that all edges incident to 'a' had changed
		boost_adjacency_iterator v, v_end;

		for (tie(v, v_end) = adjacent_vertices(a.vertex, this->skeleton); v != v_end; ++v)
			if (this->visitor) this->visitor->on_changed_edge(a,Vertex_handle(*v));
	}

	//@}


};

}

}  // namespace GUDHI

#endif /* GUDHI_SKELETON_BLOCKERS_SIMPLIFIABLE_COMPLEX_H_ */
