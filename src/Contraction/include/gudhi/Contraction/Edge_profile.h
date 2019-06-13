/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONTRACTION_EDGE_PROFILE_H_
#define CONTRACTION_EDGE_PROFILE_H_

#include <ostream>

namespace Gudhi {

namespace contraction {

template<typename GeometricSimplifiableComplex> class Edge_profile {
 public:
  typedef GeometricSimplifiableComplex Complex;
  typedef typename Complex::GT GT;

  typedef typename GeometricSimplifiableComplex::Vertex_handle Vertex_handle;
  typedef typename GeometricSimplifiableComplex::Root_vertex_handle Root_vertex_handle;


  typedef typename GeometricSimplifiableComplex::Edge_handle Edge_handle;
  typedef typename GeometricSimplifiableComplex::Graph_vertex Graph_vertex;
  typedef typename GeometricSimplifiableComplex::Graph_edge Graph_edge;
  typedef typename GeometricSimplifiableComplex::Point Point;

  Edge_profile(GeometricSimplifiableComplex& complex, Edge_handle edge) : complex_(complex), edge_handle_(edge),
      v0_(complex_.first_vertex(edge_handle_)), v1_(complex_.second_vertex(edge_handle_)) {
    assert(complex_.get_address(complex_[edge_handle_].first()));
    assert(complex_.get_address(complex_[edge_handle_].second()));
    assert(complex_.contains_edge(v0_handle(), v1_handle()));
    assert(v0_handle() != v1_handle());
  }

  virtual ~Edge_profile() { }

  GeometricSimplifiableComplex& complex() const {
    return complex_;
  }

  Edge_handle edge_handle() const {
    return edge_handle_;
  }

  Graph_edge& edge() const {
    return complex_[edge_handle_];
  }

  Graph_vertex& v0() const {
    return complex_[v0_handle()];
  }

  Graph_vertex& v1() const {
    return complex_[v1_handle()];
  }

  Vertex_handle v0_handle() const {
    return v0_;
    // Root_vertex_handle root = complex_[edge_handle_].first();
    // assert(complex_.get_address(root));
    // return *complex_.get_address(root);
  }

  Vertex_handle v1_handle() const {
    return v1_;
    // Root_vertex_handle root = complex_[edge_handle_].second();
    // assert(complex_.get_address(root));
    // return *complex_.get_address(root);
  }

  const Point& p0() const {
    return complex_.point(v0_handle());
  }

  const Point& p1() const {
    return complex_.point(v1_handle());
  }

  friend std::ostream& operator<<(std::ostream& o, const Edge_profile& v) {
    return o << "v0:" << v.v0_handle() << " v1:" << v.v1_handle();
  }

 private:
  GeometricSimplifiableComplex& complex_;

  Edge_handle edge_handle_;

  Vertex_handle v0_;

  Vertex_handle v1_;
};

template<typename EdgeProfile> class Edge_profile_factory {
 public:
  typedef typename EdgeProfile::Edge_handle Edge_handle_;
  typedef typename EdgeProfile::Complex Complex_;

  virtual EdgeProfile make_profile(
                                   Complex_& complex,
                                   Edge_handle_ edge) const {
    return EdgeProfile(complex, edge);
  }

  virtual ~Edge_profile_factory() { }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_EDGE_PROFILE_H_
