/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_SIBLINGS_H_
#define SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_SIBLINGS_H_

#include "boost/container/flat_map.hpp"
#include "Simplex_tree_node_explicit_storage.h"

#include <utility>
#include <vector>

namespace Gudhi {

/* \addtogroup simplex_tree
 * Represents a set of node of a Simplex_tree that share the same parent.
 * @{
 */

/* \brief Data structure to store a set of nodes in a SimplexTree sharing
 * the same parent node.*/
template<class SimplexTree, class MapContainer>
class Simplex_tree_siblings {
// private:
//  friend SimplexTree;
 public:
  template<class T> friend class Simplex_tree_simplex_vertex_iterator;
  template<class T> friend class Simplex_tree_boundary_simplex_iterator;
  template<class T> friend class Simplex_tree_complex_simplex_iterator;
  template<class T> friend class Simplex_tree_skeleton_simplex_iterator;

  typedef typename SimplexTree::Vertex_handle Vertex_handle;
  typedef typename SimplexTree::Filtration_value Filtration_value;
  typedef typename SimplexTree::Node Node;
  typedef MapContainer Dictionary;
  typedef typename MapContainer::iterator Dictionary_it;

  /* Default constructor.*/
  Simplex_tree_siblings()
      : oncles_(NULL),
        parent_(-1),
        members_() {
  }

  /* Constructor with values.*/
  Simplex_tree_siblings(Simplex_tree_siblings * oncles, Vertex_handle parent)
      : oncles_(oncles),
        parent_(parent),
        members_() {
  }

  /* \brief Constructor with initialized set of members.
   *
   * 'members' must be sorted and unique.*/
  template<typename RandomAccessVertexRange>
  Simplex_tree_siblings(Simplex_tree_siblings * oncles, Vertex_handle parent, const RandomAccessVertexRange & members)
      : oncles_(oncles),
        parent_(parent),
        members_(boost::container::ordered_unique_range, members.begin(),
                 members.end()) {
    for (auto& map_el : members_) {
      map_el.second.assign_children(this);
	}
  }

  /*
   * \brief Inserts a Node in the set of siblings nodes.
   *
   * If already present, assigns the minimal filtration value 
   * between input filtration_value and the value already 
   * present in the node.
   */
  void insert(Vertex_handle v, Filtration_value filtration_value) {
    auto ins = members_.emplace(v, Node(this, filtration_value));
    if (!ins.second && filtration(ins.first) > filtration_value)
      ins.first->second.assign_filtration(filtration_value);
  }

  Dictionary_it find(Vertex_handle v) {
    return members_.find(v);
  }

  Simplex_tree_siblings * oncles() {
    return oncles_;
  }

  Vertex_handle parent() const {
    return parent_;
  }

  Dictionary & members() {
    return members_;
  }

  size_t size() const {
    return members_.size();
  }

  Simplex_tree_siblings * oncles_;
  Vertex_handle parent_;
  Dictionary members_;
};

/* @} */  // end addtogroup simplex_tree
}  // namespace Gudhi

#endif  // SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_SIBLINGS_H_
