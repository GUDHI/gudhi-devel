/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SKELETON_BLOCKER_SKELETON_BLOCKER_SIMPLE_GEOMETRIC_TRAITS_H_
#define SKELETON_BLOCKER_SKELETON_BLOCKER_SIMPLE_GEOMETRIC_TRAITS_H_

#include <gudhi/Skeleton_blocker/Skeleton_blocker_simple_traits.h>

#include <string>
#include <sstream>

namespace Gudhi {

namespace skeleton_blocker {

/**
 * @extends SkeletonBlockerGeometricDS
 * @ingroup skbl
 * @brief Simple traits that is a model of SkeletonBlockerGeometricDS and
 * can be passed as a template argument to Skeleton_blocker_geometric_complex
 */
template<typename GeometryTrait>
struct Skeleton_blocker_simple_geometric_traits :
public Skeleton_blocker_simple_traits {
 public:
  typedef GeometryTrait GT;
  typedef typename GT::Point Point;
  typedef typename Skeleton_blocker_simple_traits::Root_vertex_handle Root_vertex_handle;
  typedef typename Skeleton_blocker_simple_traits::Graph_vertex Simple_vertex;

  /**
   * @brief Vertex with a point attached.
   */
  class Simple_geometric_vertex : public Simple_vertex {
    template<class ComplexGeometricTraits> friend class Skeleton_blocker_geometric_complex;
   private:
    Point point_;

    Point& point() {
      return point_;
    }

    const Point& point() const {
      return point_;
    }
  };

  class Simple_geometric_edge :
  public Skeleton_blocker_simple_traits::Graph_edge {
    int index_;
   public:
    Simple_geometric_edge()
        : Skeleton_blocker_simple_traits::Graph_edge(),
        index_(-1) { }

    /**
     * @brief Allows to modify the index of the edge.
     * The indices of the edge are used to store heap information
     * in the edge contraction algorithm.
     */
    int& index() {
      return index_;
    }

    int index() const {
      return index_;
    }
  };

  typedef Simple_geometric_vertex Graph_vertex;
  typedef Skeleton_blocker_simple_traits::Graph_edge Graph_edge;
};

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_SKELETON_BLOCKER_SIMPLE_GEOMETRIC_TRAITS_H_
