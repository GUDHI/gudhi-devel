/*
 * Skeleton_blockers_simplices_iterators.h
 *
 *  Created on: Sep 26, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_KELETON_BLOCKERS_SIMPLICES_ITERATORS_H_
#define GUDHI_SKELETON_BLOCKERS_SIMPLICES_ITERATORS_H_

#include <memory>
#include <list>
#include <iostream>
#include "gudhi/Utils.h"
#include "boost/iterator/iterator_facade.hpp"


#include "gudhi/Skeleton_blocker_link_complex.h"
#include "gudhi/Skeleton_blocker/Skeleton_blocker_link_superior.h"




namespace Gudhi {


namespace skbl {


/**
 * Link may be Skeleton_blocker_link_complex<SkeletonBlockerComplex> to iterate over all
 * simplices around a vertex OR
 * Skeleton_blocker_superior_link_complex<SkeletonBlockerComplex> to iterate over all
 * superior vertices around a vertex.
 * The iteration is done by computing a trie with the link and doing a breadth-first traversal
 * of the trie.
 */
template<typename SkeletonBlockerComplex,typename Link>
class Simplex_around_vertex_iterator :
		public boost::iterator_facade < Simplex_around_vertex_iterator<SkeletonBlockerComplex,Link>
, typename SkeletonBlockerComplex::Simplex_handle
, boost::forward_traversal_tag
, typename SkeletonBlockerComplex::Simplex_handle
>
{
	friend class boost::iterator_core_access;
	typedef SkeletonBlockerComplex Complex;
	typedef typename Complex::Vertex_handle Vertex_handle;
	typedef typename Complex::Edge_handle Edge_handle;
	typedef typename Complex::Simplex_handle Simplex_handle;


	typedef typename Link::Vertex_handle Link_vertex_handle;
	// Link_vertex_handle == Complex_Vertex_handle but this renaming helps avoiding confusion


private:

	struct Trie{
		Vertex_handle v;
		std::vector<std::shared_ptr<Trie> > childs;
	private:
		const Trie* parent_;
	public:


		//std::list<std::unique_ptr<Trie> > childs; -> use of deleted function

		Trie():parent_(0){}
		Trie(Vertex_handle v_):v(v_),parent_(0){
		}

		Trie(Vertex_handle v_,Trie* parent):v(v_),parent_(parent){
		}


		bool operator==(const Trie& other) const{
			return (v == other.v) ;
		}

		void add_child(Trie* child){
			if(child){
				std::shared_ptr<Trie> ptr_to_add(child);
				childs.push_back(ptr_to_add);
				child->parent_ = this;
			}
		}


		friend std::ostream& operator<<(std::ostream& stream, const Trie& trie){
			stream<< "T( "<< trie.v<< " ";
			for(auto t : trie.childs)
				stream << *t ;
			stream<<")";
			return stream;
		}

		// goes to the root in the trie to consitute simplex
		void add_vertices_up_to_the_root(Simplex_handle& res) const{
			res.add_vertex(v);
			if(parent_)
				parent_->add_vertices_up_to_the_root(res);
		}

		Simplex_handle simplex() const{
			Simplex_handle res;
			add_vertices_up_to_the_root(res);
			return res;
		}

		bool is_leaf() const{
			return childs.empty();
		}

		bool is_root() const{
			return parent_==0;
		}

		const Trie* parent() {
			return parent_;
		}

		void remove_leaf() {
			assert(is_leaf);
			if(!is_root())
				parent_->childs.erase(this);
		}

	private:


	public:

		Trie* go_bottom_left(){
			if(is_leaf())
				return this;
			else
				return (*childs.begin())->go_bottom_left();
		}

	};

private:
	const Complex* complex;
	Vertex_handle v;
	std::shared_ptr<Link> link_v;
	std::shared_ptr<Trie> trie;
	std::list<Trie*> nodes_to_be_seen; // todo regrouper

public:
	Simplex_around_vertex_iterator():complex(0){
	}

	Simplex_around_vertex_iterator(const Complex* complex_,Vertex_handle v_):
		complex(complex_),
		v(v_),
		link_v(new Link(*complex_,v_)),
		trie(new Trie(v_)){
		compute_trie_and_nodes_to_be_seen();
	}

	// todo avoid useless copy
	// todo currently just work if copy begin iterator
	Simplex_around_vertex_iterator(const Simplex_around_vertex_iterator& other):
		complex(other.complex),
		v(other.v),
		link_v(other.link_v),
		trie(other.trie),
		nodes_to_be_seen(other.nodes_to_be_seen)
	{
		if(!other.is_end()){
		}
	}

	/**
	 * returns an iterator to the end
	 */
	Simplex_around_vertex_iterator(const Complex* complex_,Vertex_handle v_,bool end):
		complex(complex_),
		v(v_){
		set_end();
	}

private:


	void compute_trie_and_nodes_to_be_seen(){
		// once we go through every simplices passing through v0
		// we remove v0. That way, it prevents from passing a lot of times
		// though edges leaving v0.
		// another solution would have been to provides an adjacency iterator
		// to superior vertices that avoids lower ones.
		while(!link_v->empty()){
			auto v0 = *(link_v->vertex_range().begin());
			trie->add_child(build_trie(v0,trie.get()));
			link_v->remove_vertex(v0);
		}
		nodes_to_be_seen.push_back(trie.get());
	}

	Trie* build_trie(Link_vertex_handle link_vh,Trie* parent){
		Trie* res = new Trie(parent_vertex(link_vh),parent);
		for(Link_vertex_handle nv : link_v->vertex_range(link_vh)) {
			if(link_vh < nv){
				Simplex_handle simplex_node_plus_nv(res->simplex());
				simplex_node_plus_nv.add_vertex(parent_vertex(nv));
				if(complex->contains(simplex_node_plus_nv)){
					res->add_child(build_trie(nv,res));
				}
			}
		}
		return res;
	}

	bool is_node_in_complex(Trie* trie){
		return true;
	}

	Vertex_handle parent_vertex(Link_vertex_handle link_vh) const{
		return complex->convert_handle_from_another_complex(*link_v,link_vh);
	}



public:

	friend std::ostream& operator<<(std::ostream& stream, const Simplex_around_vertex_iterator& savi){
		stream<< savi.trie<< std::endl; ;
		stream << "("<<savi.nodes_to_be_seen.size()<<") nodes to see\n";
		return stream;
	}

	bool equal(const Simplex_around_vertex_iterator& other) const{
		bool same_complex = (complex == other.complex);
		if(!same_complex)
			return false;

		if(!(v == other.v))
			return false;

		bool both_empty = nodes_to_be_seen.empty() && other.nodes_to_be_seen.empty();
		if(both_empty)
			return true;

		bool both_non_empty = !nodes_to_be_seen.empty() && !other.nodes_to_be_seen.empty();

		if(!both_non_empty) return false; //one is empty the other is not

		bool same_node = (**(nodes_to_be_seen.begin()) == **(other.nodes_to_be_seen.begin()));
		return same_node;
	}

	void increment(){
		assert(!is_end());
		Trie* first_node = nodes_to_be_seen.front();

		nodes_to_be_seen.pop_front();

		for(auto childs : first_node->childs){
			nodes_to_be_seen.push_back(childs.get());
		}

	}

	Simplex_handle dereference() const{
		assert(!nodes_to_be_seen.empty());
		Trie* first_node = nodes_to_be_seen.front();
		return first_node->simplex();
	}

private:
	void set_end(){
		nodes_to_be_seen.clear();
	}

	bool is_end() const{
		return nodes_to_be_seen.empty();
	}
};



template<typename SkeletonBlockerComplex>
class Simplex_iterator :
		public boost::iterator_facade < Simplex_iterator<SkeletonBlockerComplex>
, typename SkeletonBlockerComplex::Simplex_handle
, boost::forward_traversal_tag
, typename SkeletonBlockerComplex::Simplex_handle
>
{
	typedef Skeleton_blocker_link_superior<SkeletonBlockerComplex> Link;

	friend class boost::iterator_core_access;

	template<class SkBlDS> friend class Skeleton_blocker_complex;


	typedef SkeletonBlockerComplex Complex;
	typedef typename Complex::Vertex_handle Vertex_handle;
	typedef typename Complex::Edge_handle Edge_handle;
	typedef typename Complex::Simplex_handle Simplex_handle;

	typedef typename Complex::CVI CVI;


	typedef typename Link::Vertex_handle Link_vertex_handle;

private:

	const Complex* complex_;
	CVI current_vertex_;

	typedef Simplex_around_vertex_iterator<SkeletonBlockerComplex,Link> SAVI;
	SAVI current_simplex_around_current_vertex_;
	SAVI simplices_around_current_vertex_end_;


public:
	Simplex_iterator():complex_(0){}

	// should not be called if the complex is empty
	Simplex_iterator(const Complex* complex):
		complex_(complex),
		current_vertex_(complex->vertex_range().begin()),
		current_simplex_around_current_vertex_(complex,*current_vertex_),
		simplices_around_current_vertex_end_(complex,*current_vertex_,true)
	{
		assert(!complex->empty());
	}

private:
	// todo return to private
	Simplex_iterator(const Complex* complex,bool end):
		complex_(complex)
	{
		set_end();
	}

public:

	Simplex_iterator(const Simplex_iterator& other)
	:
		complex_(other.complex_),
		current_vertex_(other.current_vertex_),
		current_simplex_around_current_vertex_(other.current_simplex_around_current_vertex_),
		simplices_around_current_vertex_end_(other.simplices_around_current_vertex_end_)
	{
	}

	friend Simplex_iterator make_begin_iterator(const Complex* complex){
		if(complex->empty())
			return make_end_simplex_iterator(complex);
		else
			return Simplex_iterator(complex);
	}

	friend Simplex_iterator make_end_simplex_iterator(const Complex* complex){
		return Simplex_iterator(complex,true);
	}

	bool equal(const Simplex_iterator& other) const{
		if(complex_!=other.complex_) return false;
		if(current_vertex_!=other.current_vertex_) return false;
		if(is_end() && other.is_end()) return true;
		if(current_simplex_around_current_vertex_ != other.current_simplex_around_current_vertex_)
			return false;
		return true;
	}

	void increment(){
		if(current_simplex_around_current_vertex_!= simplices_around_current_vertex_end_){
			current_simplex_around_current_vertex_.increment();
			if( current_simplex_around_current_vertex_== simplices_around_current_vertex_end_)
				goto_next_vertex();
		}
		else{
			goto_next_vertex();
		}
	}

	void goto_next_vertex(){
		current_vertex_.increment();
		if(!is_end()){
			current_simplex_around_current_vertex_= SAVI(complex_,*current_vertex_);
			simplices_around_current_vertex_end_ = SAVI(complex_,*current_vertex_,true);
		}
	}

	Simplex_handle dereference() const{
		return current_simplex_around_current_vertex_.dereference();
	}

private:
	void set_end(){
		current_vertex_ = complex_->vertex_range().end();
	}

	bool is_end() const{
		return (current_vertex_ == complex_->vertex_range().end());
	}
};


}

} // namespace GUDHI





#endif /* SKELETON_BLOCKERS_SIMPLICES_ITERATORS_H_ */
