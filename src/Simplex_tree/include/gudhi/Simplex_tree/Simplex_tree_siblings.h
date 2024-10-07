/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_SIMPLEX_TREE_SIBLINGS_H_
#define SIMPLEX_TREE_SIMPLEX_TREE_SIBLINGS_H_

#include <gudhi/Simplex_tree/simplex_tree_options.h>
#include <gudhi/Simplex_tree/Simplex_tree_node_explicit_storage.h>

#include <boost/container/flat_map.hpp>

namespace Gudhi {

/** \addtogroup simplex_tree
 * Represents a set of node of a Simplex_tree that share the same parent.
 * @{
 */

/** \brief Data structure to store a set of nodes in a SimplexTree sharing
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
  template<class T> friend class Simplex_tree_boundary_opposite_vertex_simplex_iterator;

  typedef typename SimplexTree::Vertex_handle Vertex_handle;
  typedef typename SimplexTree::Filtration_value Filtration_value;
  typedef typename SimplexTree::Node Node;
  typedef MapContainer Dictionary;
  typedef typename MapContainer::iterator Dictionary_it;

  /* Default constructor.*/
  Simplex_tree_siblings()
      : oncles_(nullptr),
        parent_(-1),
        members_() {
  }

  /* Constructor with values.*/
  Simplex_tree_siblings(Simplex_tree_siblings * oncles, Vertex_handle parent)
      : oncles_(oncles),
        parent_(parent),
        members_() {
  }

  /** \brief Constructor with initialized set of members.
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

  void erase(const Dictionary_it iterator) {
    members_.erase(iterator);
  }

  Simplex_tree_siblings * oncles_;
  Vertex_handle parent_;
  Dictionary members_;
};

/** @}*/  // end addtogroup simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SIMPLEX_TREE_SIBLINGS_H_
