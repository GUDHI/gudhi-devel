/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_HOOKS_SIMPLEX_BASE_H_
#define SIMPLEX_TREE_HOOKS_SIMPLEX_BASE_H_

#include <boost/intrusive/list.hpp>

namespace Gudhi {
  /** \brief Data structure to put all simplex tree nodes with same label into a list.
   *
   * Allows one to access all subtrees of the simplex tree rooted at a node with a given label.
   * Used in particular for fast cofaces location, and fast insertion and deletion of edges in a flag complex.
   *
   * Only if SimplexTreeOptions::link_nodes_by_label is true, otherwise store nothing.
   */

  // no hook
  struct Hooks_simplex_base_dummy {};

  // todo on Hooks_simplex_base_link_nodes:
  // make the class movable but not copiable
  // DONE : translate the boost macros into C++11 syntax (boost independent)
  // update the Node concept and the doc
  struct Hooks_simplex_base_link_nodes {
   public:
    Hooks_simplex_base_link_nodes() {}
    // the copy constructor, inherited by the Node class, exchanges hooks,
    // and make the ones of this invalid.
    // this is used when stored in a map like DS, using copies when
    // performing insertions and rebalancing of the rbtree
    Hooks_simplex_base_link_nodes(const Hooks_simplex_base_link_nodes& other) {
      list_max_vertex_hook_.swap_nodes(other.list_max_vertex_hook_);
    }
    // copy assignment
    Hooks_simplex_base_link_nodes& operator=(const Hooks_simplex_base_link_nodes& other) {
      list_max_vertex_hook_.swap_nodes(other.list_max_vertex_hook_);
      return *this;
    }
    // move constructor
    Hooks_simplex_base_link_nodes(Hooks_simplex_base_link_nodes&& other) {
      list_max_vertex_hook_.swap_nodes(other.list_max_vertex_hook_);
    }
    // move assignment
    Hooks_simplex_base_link_nodes& operator=(Hooks_simplex_base_link_nodes&& other) {
      list_max_vertex_hook_.swap_nodes(other.list_max_vertex_hook_);
      return *this;
    }
    void unlink_hooks() { list_max_vertex_hook_.unlink(); }
    ~Hooks_simplex_base_link_nodes() {}  // unlink_hooks(); }

    typedef boost::intrusive::list_member_hook<  // allows .unlink()
        boost::intrusive::link_mode<boost::intrusive::auto_unlink>>
        Member_hook_t;

    mutable Member_hook_t list_max_vertex_hook_;
  };

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_HOOKS_SIMPLEX_BASE_H_
