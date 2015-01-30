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
#ifndef SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_SKELETON_BLOCKER_SKELETON_BLOCKER_SIMPLE_TRAITS_H_
#define SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_SKELETON_BLOCKER_SKELETON_BLOCKER_SIMPLE_TRAITS_H_

#include <string>
#include <sstream>
#include "Skeleton_blocker_simplex.h"

namespace Gudhi {

namespace skbl {

/**
 * @extends SkeletonBlockerDS
 * @ingroup skbl
 * @brief Simple traits that is a model of SkeletonBlockerDS and
 * can be passed as a template argument to Skeleton_blocker_complex
 */
struct Skeleton_blocker_simple_traits {
  /**
   * @brief Global and local handle similar to <a href="http://www.boost.org/doc/libs/1_38_0/libs/graph/doc/subgraph.html">boost subgraphs</a>.
   * Vertices are stored in a vector.
   * For the root simplicial complex, the local and global descriptors are the same.
   * For a subcomplex L and one of its vertices 'v', the local descriptor of 'v' is its position in
   * the vertex vector of the subcomplex L whereas its global descriptor is the position of 'v'
   * in the vertex vector of the root simplicial complex.
   */
  struct Root_vertex_handle {
    typedef int boost_vertex_handle;
    explicit Root_vertex_handle(boost_vertex_handle val = -1)
        : vertex(val) {
    }
    boost_vertex_handle vertex;

    bool operator!=(const Root_vertex_handle& other) const {
      return !(this->vertex == other.vertex);
    }

    bool operator==(const Root_vertex_handle& other) const {
      return this->vertex == other.vertex;
    }

    bool operator<(const Root_vertex_handle& other) const {
      return this->vertex < other.vertex;
    }

    friend std::ostream& operator <<(std::ostream& o,
                                     const Root_vertex_handle & v) {
      o << v.vertex;
      return o;
    }
  };

  struct Vertex_handle {
    typedef int boost_vertex_handle;
    explicit Vertex_handle(boost_vertex_handle val = -1)
        : vertex(val) {
    }

    operator int() const { return (int)vertex; }

    boost_vertex_handle vertex;

    bool operator==(const Vertex_handle& other) const {
      return this->vertex == other.vertex;
    }

    bool operator!=(const Vertex_handle& other) const {
      return this->vertex != other.vertex;
    }

    bool operator<(const Vertex_handle& other) const {
      return this->vertex < other.vertex;
    }

    friend std::ostream& operator <<(std::ostream& o, const Vertex_handle & v) {
      o << v.vertex;
      return o;
    }
  };

  class Graph_vertex {
    bool is_active_;
    Root_vertex_handle id_;

   public:
    virtual ~Graph_vertex() {
    }

    void activate() {
      is_active_ = true;
    }
    void deactivate() {
      is_active_ = false;
    }
    bool is_active() const {
      return is_active_;
    }
    void set_id(Root_vertex_handle i) {
      id_ = i;
    }
    Root_vertex_handle get_id() const {
      return id_;
    }

    virtual std::string to_string() const {
      std::ostringstream res;
      res << id_;
      return res.str();
    }

    friend std::ostream& operator <<(std::ostream& o, const Graph_vertex & v) {
      o << v.to_string();
      return o;
    }
  };

  class Graph_edge {
    Root_vertex_handle a_;
    Root_vertex_handle b_;
    int index_;

   public:
    Graph_edge()
        : a_(-1),
          b_(-1),
          index_(-1) {
    }

    int& index() {
      return index_;
    }
    int index() const {
      return index_;
    }

    void setId(Root_vertex_handle a, Root_vertex_handle b) {
      a_ = a;
      b_ = b;
    }

    Root_vertex_handle first() const {
      return a_;
    }

    Root_vertex_handle second() const {
      return b_;
    }

    friend std::ostream& operator <<(std::ostream& o, const Graph_edge & v) {
      o << "(" << v.a_ << "," << v.b_ << " - id = " << v.index();
      return o;
    }
  };
};

}  // namespace skbl

}  // namespace Gudhi

#endif  // SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_SKELETON_BLOCKER_SKELETON_BLOCKER_SIMPLE_TRAITS_H_
