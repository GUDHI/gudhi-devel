/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
// TODO(DS): the index is just for contraction, to remove

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
