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

#ifndef SKELETON_BLOCKER_ITERATORS_SKELETON_BLOCKERS_EDGES_ITERATORS_H_
#define SKELETON_BLOCKER_ITERATORS_SKELETON_BLOCKERS_EDGES_ITERATORS_H_

#include <boost/iterator/iterator_facade.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <utility>  // for pair<>

namespace Gudhi {

namespace skeleton_blocker {

template<typename SkeletonBlockerComplex>
class Edge_around_vertex_iterator : public boost::iterator_facade <Edge_around_vertex_iterator<SkeletonBlockerComplex>
, typename SkeletonBlockerComplex::Edge_handle const, boost::forward_traversal_tag
, typename SkeletonBlockerComplex::Edge_handle const> {
  friend class boost::iterator_core_access;

  typedef SkeletonBlockerComplex Complex;
  typedef typename Complex::boost_adjacency_iterator boost_adjacency_iterator;
  typedef typename Complex::Vertex_handle Vertex_handle;
  typedef typename Complex::Edge_handle Edge_handle;

 private:
  const Complex* complex;
  Vertex_handle v;

  boost_adjacency_iterator current_;
  boost_adjacency_iterator end_;

 public:
  Edge_around_vertex_iterator() : complex(NULL) { }

  Edge_around_vertex_iterator(const Complex* complex_, Vertex_handle v_) :
      complex(complex_),
      v(v_) {
    tie(current_, end_) = adjacent_vertices(v.vertex, complex->skeleton);
  }

  /**
   * returns an iterator to the end
   */
  Edge_around_vertex_iterator(const Complex* complex_, Vertex_handle v_, int end) :
      complex(complex_),
      v(v_) {
    tie(current_, end_) = adjacent_vertices(v.vertex, complex->skeleton);
    set_end();
  }

  bool equal(const Edge_around_vertex_iterator& other) const {
    return (complex == other.complex) && (v == other.v) && (current_ == other.current_) && (end_ == other.end_);
  }

  void increment() {
    if (current_ != end_)
      ++current_;
  }

  Edge_handle dereference() const {
    return *(*complex)[std::make_pair(v, static_cast<Vertex_handle> (*current_))];
  }

 private:
  // remove this ugly hack
  void set_end() {
    current_ = end_;
  }
};

/**
 *@brief Iterator on the edges of a simplicial complex.
 *
 */
template<typename SkeletonBlockerComplex>
class Edge_iterator : public boost::iterator_facade <Edge_iterator<SkeletonBlockerComplex>
, typename SkeletonBlockerComplex::Edge_handle const
, boost::forward_traversal_tag
, typename SkeletonBlockerComplex::Edge_handle const> {
  friend class boost::iterator_core_access;

 public:
  typedef SkeletonBlockerComplex Complex;
  typedef typename Complex::boost_edge_iterator boost_edge_iterator;
  typedef typename Complex::Edge_handle Edge_handle;

  const Complex* complex;
  std::pair<boost_edge_iterator, boost_edge_iterator> edge_iterator;

  Edge_iterator() : complex(NULL) { }

  Edge_iterator(const SkeletonBlockerComplex* complex_) :
      complex(complex_),
      edge_iterator(boost::edges(complex_->skeleton)) { }

  /**
   * return an iterator to the end
   */
  Edge_iterator(const SkeletonBlockerComplex* complex_, int end) :
      complex(complex_),
      edge_iterator(boost::edges(complex_->skeleton)) {
    edge_iterator.first = edge_iterator.second;
  }

  bool equal(const Edge_iterator& other) const {
    return (complex == other.complex) && (edge_iterator == other.edge_iterator);
  }

  void increment() {
    if (edge_iterator.first != edge_iterator.second) {
      ++(edge_iterator.first);
    }
  }

  Edge_handle dereference() const {
    return (*(edge_iterator.first));
  }
};

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_ITERATORS_SKELETON_BLOCKERS_EDGES_ITERATORS_H_
