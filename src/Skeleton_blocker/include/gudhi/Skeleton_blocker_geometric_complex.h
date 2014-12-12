/*
 * Skeleton_blocker_geometric_complex.h
 *
 *  Created on: Feb 11, 2014
 *      Author: dsalinas
 */

#ifndef SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_
#define SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_


#include "gudhi/Utils.h"
#include "gudhi/Skeleton_blocker_simplifiable_complex.h"
#include "gudhi/Skeleton_blocker/Skeleton_blocker_sub_complex.h"



namespace Gudhi{

namespace skbl {

/**
 * @brief Class that represents a geometric complex that can be simplified.
 * The class allows access to points of vertices.
 *
 */
template<typename SkeletonBlockerGeometricDS>
class Skeleton_blocker_geometric_complex : public Skeleton_blocker_simplifiable_complex<SkeletonBlockerGeometricDS>
{
public:
	typedef typename SkeletonBlockerGeometricDS::GT GT;

	typedef Skeleton_blocker_simplifiable_complex<SkeletonBlockerGeometricDS> SimplifiableSkeletonblocker;

	typedef typename SimplifiableSkeletonblocker::Vertex_handle Vertex_handle;
	typedef typename SimplifiableSkeletonblocker::Root_vertex_handle Root_vertex_handle;
	typedef typename SimplifiableSkeletonblocker::Edge_handle Edge_handle;
	typedef typename SimplifiableSkeletonblocker::Simplex_handle Simplex_handle;


	typedef typename SimplifiableSkeletonblocker::Graph_vertex Graph_vertex;

	typedef typename SkeletonBlockerGeometricDS::Point Point;


	/**
	 * @brief Add a vertex to the complex with its associated point.
	 */
	Vertex_handle add_vertex(const Point& point){
		Vertex_handle ad = SimplifiableSkeletonblocker::add_vertex();
		(*this)[ad].point() = point;
		return ad;
	}


	/**
	 * @brief Returns the Point associated to the vertex v.
	 */
	const Point& point(Vertex_handle v) const{
		assert(this->contains_vertex(v));
		return (*this)[v].point();
	}

	/**
	 * @brief Returns the Point associated to the vertex v.
	 */
	Point& point(Vertex_handle v) {
		assert(this->contains_vertex(v));
		return (*this)[v].point();
	}

	const Point& point(Root_vertex_handle global_v) const{
		Vertex_handle local_v ( (*this)[global_v]) ;
		assert(this->contains_vertex(local_v));
		return (*this)[local_v].point();
	}

	Point& point(Root_vertex_handle global_v) {
		Vertex_handle local_v ( (*this)[global_v]) ;
		assert(this->contains_vertex(local_v));
		return (*this)[local_v].point();
	}


	typedef Skeleton_blocker_link_complex<Skeleton_blocker_geometric_complex> Geometric_link;

	/**
	 * Constructs the link of 'simplex' with points coordinates.
	 */
	Geometric_link link(const Simplex_handle& simplex) const{
		Geometric_link link(*this,simplex);
		//we now add the point info
		add_points_to_link(link);
		return link;
	}

	/**
	 * Constructs the link of 'simplex' with points coordinates.
	 */
	Geometric_link link(Edge_handle edge) const{
		Geometric_link link(*this,edge);
		//we now add the point info
		add_points_to_link(link);
		return link;
	}
private:
	void add_points_to_link(Geometric_link& link) const{
		for(Vertex_handle v : link.vertex_range()){
			Root_vertex_handle v_root(link.get_id(v));
			link.point(v) = (*this).point(v_root);
		}
	}
};

}

}  // namespace GUDHI

#endif /* SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_ */
