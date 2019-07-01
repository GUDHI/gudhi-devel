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

#ifndef SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_
#define SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_

#include <gudhi/Skeleton_blocker_complex.h>
#include <gudhi/Skeleton_blocker/Skeleton_blocker_sub_complex.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {

namespace skeleton_blocker {

/**
 * @brief Class that represents a geometric complex that can be simplified.
 * The class allows access to points of vertices.
 * @ingroup skbl
 */
template<typename SkeletonBlockerGeometricDS>
class Skeleton_blocker_geometric_complex :
public Skeleton_blocker_complex<SkeletonBlockerGeometricDS> {
 public:
  typedef typename SkeletonBlockerGeometricDS::GT GT;

  typedef Skeleton_blocker_complex<SkeletonBlockerGeometricDS> SimplifiableSkeletonblocker;

  typedef typename SimplifiableSkeletonblocker::Vertex_handle Vertex_handle;
  typedef typename SimplifiableSkeletonblocker::Root_vertex_handle Root_vertex_handle;
  typedef typename SimplifiableSkeletonblocker::Edge_handle Edge_handle;
  typedef typename SimplifiableSkeletonblocker::Simplex Simplex;

  typedef typename SimplifiableSkeletonblocker::Graph_vertex Graph_vertex;

  typedef typename SkeletonBlockerGeometricDS::Point Point;

  Skeleton_blocker_geometric_complex() { }

  /**
   * constructor given a list of points
   */
  template<typename PointIterator>
  explicit Skeleton_blocker_geometric_complex(int num_vertices, PointIterator begin, PointIterator end) {
    for (auto point = begin; point != end; ++point)
      add_vertex(*point);
  }

  /**
   * @brief Constructor with a list of simplices.
   * @details is_flag_complex indicates if the complex is a flag complex or not (to know if blockers have to be
   * computed or not).
   */
  template<typename SimpleHandleOutputIterator, typename PointIterator>
  Skeleton_blocker_geometric_complex(
                                     SimpleHandleOutputIterator simplex_begin, SimpleHandleOutputIterator simplex_end,
                                     PointIterator points_begin, PointIterator points_end,
                                     bool is_flag_complex = false)
      : Skeleton_blocker_complex<SkeletonBlockerGeometricDS>(simplex_begin, simplex_end, is_flag_complex) {
    unsigned current = 0;
    for (auto point = points_begin; point != points_end; ++point)
      (*this)[Vertex_handle(current++)].point() = Point(point->begin(), point->end());
  }

  /**
   * @brief Constructor with a list of simplices.
   * Points of every vertex are the point constructed with default constructor.
   * @details is_flag_complex indicates if the complex is a flag complex or not (to know if blockers have to be computed or not).
   */
  template<typename SimpleHandleOutputIterator>
  Skeleton_blocker_geometric_complex(
                                     SimpleHandleOutputIterator simplex_begin, SimpleHandleOutputIterator simplex_end,
                                     bool is_flag_complex = false)
      : Skeleton_blocker_complex<SkeletonBlockerGeometricDS>(simplex_begin, simplex_end, is_flag_complex) { }

  /**
   * @brief Add a vertex to the complex with a default constructed associated point.
   */
  Vertex_handle add_vertex() {
    return SimplifiableSkeletonblocker::add_vertex();
  }

  /**
   * @brief Add a vertex to the complex with its associated point.
   */
  Vertex_handle add_vertex(const Point& point) {
    Vertex_handle ad = SimplifiableSkeletonblocker::add_vertex();
    (*this)[ad].point() = point;
    return ad;
  }

  /**
   * @brief Returns the Point associated to the vertex v.
   */
  const Point& point(Vertex_handle v) const {
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

  const Point& point(Root_vertex_handle global_v) const {
    Vertex_handle local_v((*this)[global_v]);
    assert(this->contains_vertex(local_v));
    return (*this)[local_v].point();
  }

  Point& point(Root_vertex_handle global_v) {
    Vertex_handle local_v((*this)[global_v]);
    assert(this->contains_vertex(local_v));
    return (*this)[local_v].point();
  }

  typedef Skeleton_blocker_link_complex<Skeleton_blocker_geometric_complex> Geometric_link;

  /**
   * Constructs the link of 'simplex' with points coordinates.
   */
  Geometric_link link(Vertex_handle v) const {
    Geometric_link link(*this, Simplex(v));
    // we now add the point info
    add_points_to_link(link);
    return link;
  }

  /**
   * Constructs the link of 'simplex' with points coordinates.
   */
  Geometric_link link(const Simplex& simplex) const {
    Geometric_link link(*this, simplex);
    // we now add the point info
    add_points_to_link(link);
    return link;
  }

  /**
   * Constructs the link of 'simplex' with points coordinates.
   */
  Geometric_link link(Edge_handle edge) const {
    Geometric_link link(*this, edge);
    // we now add the point info
    add_points_to_link(link);
    return link;
  }

  typedef Skeleton_blocker_link_complex<Skeleton_blocker_complex<SkeletonBlockerGeometricDS>> Abstract_link;

  /**
   * Constructs the abstract link of v (without points coordinates).
   */
  Abstract_link abstract_link(Vertex_handle v) const {
    return Abstract_link(*this, Simplex(v));
  }

  /**
   * Constructs the link of 'simplex' with points coordinates.
   */
  Abstract_link abstract_link(const Simplex& simplex) const {
    return Abstract_link(*this, simplex);
  }

  /**
   * Constructs the link of 'simplex' with points coordinates.
   */
  Abstract_link abstract_link(Edge_handle edge) const {
    return Abstract_link(*this, edge);
  }

 private:
  void add_points_to_link(Geometric_link& link) const {
    for (Vertex_handle v : link.vertex_range()) {
      Root_vertex_handle v_root(link.get_id(v));
      link.point(v) = (*this).point(v_root);
    }
  }
};


template<typename SkeletonBlockerGeometricComplex, typename SimpleHandleOutputIterator, typename PointIterator>
SkeletonBlockerGeometricComplex make_complex_from_top_faces(
        SimpleHandleOutputIterator simplex_begin,
        SimpleHandleOutputIterator simplex_end,
        PointIterator points_begin,
        PointIterator points_end,
        bool is_flag_complex = false) {
  typedef SkeletonBlockerGeometricComplex SBGC;
  SkeletonBlockerGeometricComplex complex;
  unsigned current = 0;
  complex =
    make_complex_from_top_faces<SBGC>(simplex_begin, simplex_end, is_flag_complex);
  for (auto point = points_begin; point != points_end; ++point)
    // complex.point(Vertex_handle(current++)) = Point(point->begin(),point->end());
    complex.point(typename SBGC::Vertex_handle(current++)) = typename SBGC::Point(*point);
  return complex;
}

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_
