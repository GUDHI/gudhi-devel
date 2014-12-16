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
#ifndef GUDHI_SKELETON_BLOCKERS_SIMPLE_GEOMETRIC_TRAITS_H_
#define GUDHI_SKELETON_BLOCKERS_SIMPLE_GEOMETRIC_TRAITS_H_

#include <string>
#include <sstream>
#include "gudhi/Skeleton_blocker/Skeleton_blocker_simple_traits.h"

namespace Gudhi{

namespace skbl{


/**
 * @extends SkeletonBlockerGeometricDS
 * @ingroup skbl
 * @brief Simple traits that is a model of SkeletonBlockerGeometricDS and
 * can be passed as a template argument to Skeleton_blocker_geometric_complex
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
