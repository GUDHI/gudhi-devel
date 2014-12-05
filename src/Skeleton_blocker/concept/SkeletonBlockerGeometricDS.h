/*
 * SkeletonBlockerGeometricDS.h
 *
 *  Created on: Feb 20, 2014 
 *      Author: David Salinas
 *  Copyright 2013 INRIA. All rights reserved
 */

#ifndef GUDHI_SKELETONBLOCKERGEOMETRICDS_H_
#define GUDHI_SKELETONBLOCKERGEOMETRICDS_H_

/** \brief Concept that must be passed to
 * the template class Skeleton_blocker_geometric_complex
 *
 */
template<typename GeometryTrait>
struct SkeletonBlockerGeometricDS : public SkeletonBlockerDS
{

	/**
	 * Geometry information.
	 */
	typedef GeometryTrait GT ;

	/**
	 * Type of point (should be the same as GT::Point).
	 */
	typedef typename GeometryTrait::Point Point;

	/**
	 * @brief Vertex that stores a point.
	 */
	class Graph_vertex : public SkeletonBlockerDS::Graph_vertex{
	public:
		/**
		 * @brief Access to the point.
		 */
		Point& point(){	return point_; }
		/**
		 * @brief Access to the point.
		 */
		const Point& point() const {	return point_; }
	};

	/**
	 * @brief Edge that allows to access to an index.
	 * The indices of the edges are used to store heap information
	 * in the edge contraction algorithm.
	 */
	class Graph_Edge : public SkeletonBlockerDS::Graph_edge{
	public:
		/**
		 * @brief Access to the index.
		 */
		int& index();
		/**
		 * @brief Access to the index.
		 */
		int index();
	};
};



#endif /* GUDHI_SKELETONBLOCKERGEOMETRICDS_H_ */
