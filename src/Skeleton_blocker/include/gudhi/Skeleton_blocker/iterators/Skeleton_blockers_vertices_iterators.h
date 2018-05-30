/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
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

#ifndef SKELETON_BLOCKER_ITERATORS_SKELETON_BLOCKERS_VERTICES_ITERATORS_H_
#define SKELETON_BLOCKER_ITERATORS_SKELETON_BLOCKERS_VERTICES_ITERATORS_H_

#include <boost/iterator/iterator_facade.hpp>

#include <utility>  // for pair<>

namespace Gudhi {

namespace skeleton_blocker {

/**
 *@brief Iterator on the vertices of a simplicial complex
 *@remark The increment operator go to the next active vertex.
 *@remark Incrementation increases Vertex_handle.
 */
template<typename SkeletonBlockerComplex>
class Vertex_iterator : public boost::iterator_facade< Vertex_iterator <SkeletonBlockerComplex>
, typename SkeletonBlockerComplex::Vertex_handle const
, boost::forward_traversal_tag
, typename SkeletonBlockerComplex::Vertex_handle const> {
  friend class boost::iterator_core_access;

  typedef typename SkeletonBlockerComplex::boost_vertex_iterator boost_vertex_iterator;
  typedef typename SkeletonBlockerComplex::Vertex_handle Vertex_handle;
 private:
  const SkeletonBlockerComplex* complex;
  std::pair<boost_vertex_iterator, boost_vertex_iterator> vertexIterator;


 public:
  Vertex_iterator() : complex(NULL) { }

  Vertex_iterator(const SkeletonBlockerComplex* complex_) :
      complex(complex_),
      vertexIterator(vertices(complex_->skeleton)) {
    if (!finished() && !is_active()) {
      goto_next_valid();
    }
  }

  /**
   * return an iterator to the end.
   */
  Vertex_iterator(const SkeletonBlockerComplex* complex_, int end) :
      complex(complex_), vertexIterator(vertices(complex_->skeleton)) {
    vertexIterator.first = vertexIterator.second;
  }

 public:
  void increment() {
    goto_next_valid();
  }

  Vertex_handle dereference() const {
    return (Vertex_handle(*(vertexIterator.first)));
  }

  bool equal(const Vertex_iterator& other) const {
    return vertexIterator == other.vertexIterator && complex == other.complex;
  }

  bool operator<(const Vertex_iterator& other) const {
    return dereference() < other.dereference();
  }

 private:
  bool finished() const {
    return vertexIterator.first == vertexIterator.second;
  }

  void goto_next_valid() {
    ++vertexIterator.first;
    if (!finished() && !is_active()) {
      goto_next_valid();
    }
  }

  bool is_active() const {
    return ((*complex)[Vertex_handle(*vertexIterator.first)]).is_active();
  }
};

template<typename SkeletonBlockerComplex>
class Neighbors_vertices_iterator : public boost::iterator_facade < Neighbors_vertices_iterator<SkeletonBlockerComplex>
, typename SkeletonBlockerComplex::Vertex_handle const
, boost::forward_traversal_tag
, typename SkeletonBlockerComplex::Vertex_handle const> {
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
  Neighbors_vertices_iterator() : complex(NULL) { }

  Neighbors_vertices_iterator(const Complex* complex_, Vertex_handle v_) :
      complex(complex_),
      v(v_) {
    tie(current_, end_) = adjacent_vertices(v.vertex, complex->skeleton);
  }

  /**
   * returns an iterator to the end
   */
  Neighbors_vertices_iterator(const Complex* complex_, Vertex_handle v_, int end) :
      complex(complex_),
      v(v_) {
    tie(current_, end_) = adjacent_vertices(v.vertex, complex->skeleton);
    set_end();
  }

  void increment() {
    if (current_ != end_)
      ++current_;
  }

  Vertex_handle dereference() const {
    return (Vertex_handle(*current_));
  }

  bool equal(const Neighbors_vertices_iterator& other) const {
    return (complex == other.complex) && (v == other.v) && (current_ == other.current_) && (end_ == other.end_);
  }

 private:
  // TODO(DS): remove this ugly hack
  void set_end() {
    current_ = end_;
  }
};

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_ITERATORS_SKELETON_BLOCKERS_VERTICES_ITERATORS_H_
