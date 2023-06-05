/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_NODES_BY_LABEL_H_
#define SIMPLEX_TREE_NODES_BY_LABEL_H_

#include <boost/intrusive/list.hpp>

#include <unordered_map>

namespace Gudhi {

  // trivial data structure, in case the option is disallowed.
  template <typename SimplexTree>
  struct Nodes_by_label_dummy {};

  // use the Node hooks Hooks_simplex_base to put all Nodes with same label u in an
  // intrusive list.
  template <typename SimplexTree>
  struct Nodes_by_label_intrusive_list {
    Nodes_by_label_intrusive_list(){};
    // insert a Node in the hook list corresponding to its label
    void insert(typename SimplexTree::Simplex_handle sh) {
      auto it = nodes_label_to_list_.find(sh->first);
      if (it == nodes_label_to_list_.end()) {  // create a new list
        it = (nodes_label_to_list_.emplace(sh->first, new typename SimplexTree::List_max_vertex())).first;
      }
      it->second->push_back(sh->second);  // insert at the end of the list
    }

    ~Nodes_by_label_intrusive_list() {
      for (auto u_list_ptr : nodes_label_to_list_) {
        delete u_list_ptr.second;
      }
    }
    typename SimplexTree::List_max_vertex* find(typename SimplexTree::Vertex_handle v) {
      auto it_v = nodes_label_to_list_.find(v);
      if (it_v != nodes_label_to_list_.end()) {
        return it_v->second;
      } else {
        return nullptr;
      }
    }
    // unordered_map Vertex_handle v -> pointer to list of all Nodes with label v.
    std::unordered_map<typename SimplexTree::Vertex_handle, typename SimplexTree::List_max_vertex*> nodes_label_to_list_;
  };

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_NODES_BY_LABEL_H_
