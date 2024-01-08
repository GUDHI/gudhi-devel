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
#include <stdexcept>

namespace Gudhi {
  /** \brief No hook when SimplexTreeOptions::link_nodes_by_label is false.
   */
  struct Hooks_simplex_base_dummy {};

  /** \brief Data structure to put all simplex tree nodes with same label into a list.
   *
   * Allows one to access all subtrees of the simplex tree rooted at a node with a given label.
   * Used in particular for fast cofaces location, and fast insertion and deletion of edges in a flag complex.
   *
   * Only if SimplexTreeOptions::link_nodes_by_label is true.
   */
  struct Hooks_simplex_base_link_nodes {
   public:
    Hooks_simplex_base_link_nodes() {}
    // The following 2 functions are conceptually wrong, but when copying a Simplex_tree it is more convenient to have
    // them (and fix things asap after they are called).
    // copy constructor, inherited by the Node class, requires to copy hooks
    Hooks_simplex_base_link_nodes(const Hooks_simplex_base_link_nodes&)
    : list_max_vertex_hook_() { }
    // copy assignment
    Hooks_simplex_base_link_nodes& operator=(const Hooks_simplex_base_link_nodes&) {
      throw std::logic_error("Should not happen");
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

    Member_hook_t list_max_vertex_hook_;
  };

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_HOOKS_SIMPLEX_BASE_H_
