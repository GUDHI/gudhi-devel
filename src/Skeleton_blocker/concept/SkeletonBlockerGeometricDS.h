/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA
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

#ifndef CONCEPT_SKELETON_BLOCKER_SKELETONBLOCKERGEOMETRICDS_H_
#define CONCEPT_SKELETON_BLOCKER_SKELETONBLOCKERGEOMETRICDS_H_

namespace Gudhi {
namespace skeleton_blocker {

/** 
 * \brief Concept for template class of  Skeleton_blocker_geometric_complex .
 * It must specify a GeometryTrait which contains a Point definition.
 * 
 * Graph_vertex must specify how to access to a point.
 * Graph_edge must specify how to access to an index.
 * 
 */
//todo the index is just for contraction, to remove

template<typename GeometryTrait>
struct SkeletonBlockerGeometricDS : public SkeletonBlockerDS {
  /**
   * Geometry information.
   */
  typedef GeometryTrait GT;

  /**
   * Type of point (should be the same as GT::Point).
   */
  typedef typename GeometryTrait::Point Point;

  /**
   * @brief Vertex that stores a point.
   */
  class Graph_vertex : public SkeletonBlockerDS::Graph_vertex {
   public:
    /**
     * @brief Access to the point.
     */
    Point& point();
    /**
     * @brief Access to the point.
     */
    const Point& point();
  };

  /**
   * @brief Edge that allows to access to an index.
   * The indices of the edges are used to store heap information
   * in the edge contraction algorithm.
   */
  class Graph_Edge : public SkeletonBlockerDS::Graph_edge {
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

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif  // CONCEPT_SKELETON_BLOCKER_SKELETONBLOCKERGEOMETRICDS_H_
