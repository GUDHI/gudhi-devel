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
