/*
 * Skeleton_blockers_triangles_iterators.h
 *
 *  Created on: Aug 25, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_SKELETON_BLOCKERS_TRIANGLES_ITERATORS_H_
#define GUDHI_SKELETON_BLOCKERS_TRIANGLES_ITERATORS_H_

#include "boost/iterator/iterator_facade.hpp"
#include <memory>

namespace Gudhi{

namespace skbl {

//////////////////////////////////////////////////////////////////////
/**
 * \brief Iterator over the triangles that are
 * adjacent to a vertex of the simplicial complex.
 * \remark Will be removed soon -> dont look
 */
template<typename Complex,typename LinkType>
class Triangle_around_vertex_iterator : public boost::iterator_facade
< Triangle_around_vertex_iterator <Complex,LinkType>
, typename Complex::Simplex_handle const
, boost::forward_traversal_tag
, typename Complex::Simplex_handle const
>
{
	friend class boost::iterator_core_access;
	template<typename T> friend class Triangle_iterator;
private:
	typedef typename LinkType::Vertex_handle Vertex_handle;
	typedef typename LinkType::Root_vertex_handle Root_vertex_handle;
	typedef typename LinkType::Simplex_handle Simplex_handle;
	typedef Complex_edge_iterator<Complex> Complex_edge_iterator_;

	const Complex* complex_;
	Vertex_handle v_;
	std::shared_ptr<LinkType> link_;
	Complex_edge_iterator_ current_edge_;
	bool is_end_;
public:
	Triangle_around_vertex_iterator(const Complex* complex,Vertex_handle v):
		complex_(complex),v_(v),link_(new LinkType(*complex,v_)),
		current_edge_(link_->edge_range().begin()),
		is_end_(current_edge_ == link_->edge_range().end()){
	}

	/**
	 * @brief ugly hack to get an iterator to the end
	 */
	Triangle_around_vertex_iterator(const Complex* complex,Vertex_handle v,bool is_end):
		complex_(complex),v_(v),link_(0),is_end_(true){
	}

	/**
	 * @brief ugly hack to get an iterator to the end
	 */
	Triangle_around_vertex_iterator():
		complex_(0),v_(-1),link_(0),is_end_(true){
	}


	Triangle_around_vertex_iterator(const Triangle_around_vertex_iterator& other){
		v_ = other.v_;
		complex_ = other.complex_;
		is_end_ = other.is_end_;

		if(!is_end_){
			link_ = other.link_;
			current_edge_= other.current_edge_;
		}
	}

	bool equal(const Triangle_around_vertex_iterator& other) const{
		return (complex_==other.complex_) && ((finished() &&other.finished()) || current_edge_ == other.current_edge_);
	}

	Simplex_handle dereference() const{
		Root_vertex_handle v1 = (*link_)[*current_edge_].first();
		Root_vertex_handle v2 = (*link_)[*current_edge_].second();
		return Simplex_handle(v_,*(complex_->get_address(v1)),*(complex_->get_address(v2)));
	}

	void increment(){
		++current_edge_;
	}

private:
	bool finished() const{
		return is_end_ || (current_edge_ == link_->edge_range().end());
	}

};



/**
 * \brief Iterator over the triangles of the
 * simplicial complex.
 * \remark Will be removed soon -> dont look
 *
 */
template<typename SkeletonBlockerComplex>
class Triangle_iterator :
		public boost::iterator_facade<
		Triangle_iterator <SkeletonBlockerComplex>,
		typename SkeletonBlockerComplex::Simplex_handle const
		, boost::forward_traversal_tag
		, typename SkeletonBlockerComplex::Simplex_handle const
		>
{
	friend class boost::iterator_core_access;
private:
	typedef typename SkeletonBlockerComplex::Vertex_handle Vertex_handle;
	typedef typename SkeletonBlockerComplex::Root_vertex_handle Root_vertex_handle;
	typedef typename SkeletonBlockerComplex::Simplex_handle Simplex_handle;
	typedef typename SkeletonBlockerComplex::Superior_triangle_around_vertex_iterator STAVI;

	const SkeletonBlockerComplex* complex_;
	Complex_vertex_iterator<SkeletonBlockerComplex> current_vertex_;
	STAVI current_triangle_;
	bool is_end_;
public:

	/*
	 * @remark  assume that the complex is non-empty
	 */
	Triangle_iterator(const SkeletonBlockerComplex* complex):
		complex_(complex),
		current_vertex_(complex->vertex_range().begin()),
		current_triangle_(complex,*current_vertex_), // xxx this line is problematic is the complex is empty
		is_end_(false){

		assert(!complex->empty());
		gotoFirstTriangle();
	}

private:
	//goto to the first triangle or to the end if none
	void gotoFirstTriangle(){
		if(!is_finished() && current_triangle_.finished()){
			goto_next_vertex();
		}
	}
public:

	/**
	 * @brief ugly hack to get an iterator to the end
	 * @remark  assume that the complex is non-empty
	 */
	Triangle_iterator(const SkeletonBlockerComplex* complex,bool is_end):
		complex_(complex),
		current_vertex_(complex->vertex_range().end()),
		current_triangle_(), // xxx this line is problematic is the complex is empty
		is_end_(true){
	}


	Triangle_iterator& operator=(const Triangle_iterator & other){
		complex_ = other.complex_;
		Complex_vertex_iterator<SkeletonBlockerComplex> current_vertex_;
		STAVI current_triangle_;
		return *this;
	}


	bool equal(const Triangle_iterator& other) const{
		bool both_are_finished = is_finished() && other.is_finished();
		bool both_arent_finished = !is_finished() && !other.is_finished();
		// if the two iterators are not finished, they must have the same state
		return (complex_==other.complex_) &&
				(both_are_finished	||
						( (both_arent_finished) && current_vertex_ == other.current_vertex_ && current_triangle_ == other.current_triangle_));

	}

	Simplex_handle dereference() const{
		return *current_triangle_;
	}

private:

	// goto the next vertex that has a triangle pending or the
	// end vertex iterator if none exists
	void goto_next_vertex(){
		assert(current_triangle_.finished()); //we mush have consume all triangles passing through the vertex
		assert(!is_finished()); // we must not be done

		++current_vertex_;

		if(!is_finished()){
			current_triangle_ = STAVI(complex_, *current_vertex_);
			if(current_triangle_.finished())
				goto_next_vertex();
		}
	}
public:
	void increment(){
		if(!current_triangle_.finished()){
			++current_triangle_; // problem here
			if(current_triangle_.finished())
				goto_next_vertex();
		}
		else{
			assert(!is_finished());
			goto_next_vertex();
		}
	}

private:
	bool is_finished() const{
		return is_end_ || current_vertex_ == complex_->vertex_range().end();
	}
};

}

} // namespace GUDHI

#endif /* GUDHI_SKELETON_BLOCKERS_TRIANGLES_ITERATORS_H_ */
