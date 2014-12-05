/*
 * Skeleton_blocker_simple_geometric_traits.h
 *
 *  Created on: Feb 11, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_SKELETON_BLOCKERS_SIMPLE_GEOMETRIC_TRAITS_H_
#define GUDHI_SKELETON_BLOCKERS_SIMPLE_GEOMETRIC_TRAITS_H_

#include <string>
#include <sstream>
#include "gudhi/Skeleton_blocker/Skeleton_blocker_simple_traits.h"

namespace Gudhi{

namespace skbl{

/**
 * @extends SkeletonBlockerGeometricDS
 */
template<typename GeometryTrait>
struct Skeleton_blocker_simple_geometric_traits : public skbl::Skeleton_blocker_simple_traits {
public:

	typedef GeometryTrait GT;
	typedef typename GT::Point Point;
	typedef typename Skeleton_blocker_simple_traits::Root_vertex_handle Root_vertex_handle;
	typedef typename Skeleton_blocker_simple_traits::Graph_vertex Simple_vertex;

	/**
	 * @brief Vertex with a point attached.
	 */
	class Simple_geometric_vertex : public Simple_vertex{
		template<class ComplexGeometricTraits> friend class Skeleton_blocker_geometric_complex;
	private:
		Point point_;
		Point& point(){	return point_; }
		const Point& point() const {	return point_; }
	};


	class Simple_geometric_edge : public Skeleton_blocker_simple_traits::Graph_edge{
		int index_;
	public:
		Simple_geometric_edge():Skeleton_blocker_simple_traits::Graph_edge(),index_(-1){}
		/**
		 * @brief Allows to modify the index of the edge.
		 * The indices of the edge are used to store heap information
		 * in the edge contraction algorithm.
		 */
		int& index(){return index_;}
		int index() const {return index_;}
	};


	typedef Simple_geometric_vertex Graph_vertex;
	typedef Skeleton_blocker_simple_traits::Graph_edge Graph_edge;
};


}

}  // namespace GUDHI

#endif /* GUDHI_SKELETON_BLOCKERS_SIMPLE_GEOMETRIC_TRAITS_H_ */
