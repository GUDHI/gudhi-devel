/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_
#define SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_

// Empty base optimization for MSVC - https://learn.microsoft.com/en-us/cpp/cpp/empty-bases
#if _MSC_VER
 #define GUDHI_EMPTY_BASE_CLASS_OPTIMIZATION __declspec(empty_bases)
#else
 #define GUDHI_EMPTY_BASE_CLASS_OPTIMIZATION
#endif

#include <boost/core/empty_value.hpp>

namespace Gudhi {

/** \addtogroup simplex_tree
 * Represents a node of a Simplex_tree.
 * @{
 */

/** \brief Node of a simplex tree with filtration value
 * and simplex key.
 *
 * It stores explicitly its own filtration value and its own Simplex_key.
 */
template <class SimplexTree>
struct GUDHI_EMPTY_BASE_CLASS_OPTIMIZATION Simplex_tree_node_explicit_storage
    : SimplexTree::Filtration_simplex_base,
      SimplexTree::Key_simplex_base,
      SimplexTree::Hooks_simplex_base,
      boost::empty_value<typename SimplexTree::Simplex_data> {
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Filtration_value Filtration_value;
  typedef typename SimplexTree::Simplex_key Simplex_key;
  typedef typename SimplexTree::Simplex_data Simplex_data;

  // Simplex_tree_node_explicit_storage() : children_(nullptr) {}

  // Simplex_tree_node_explicit_storage(Siblings* sib, const Filtration_value& filtration) : children_(sib)
  // {
  //   this->assign_filtration(filtration);
  // }
  Simplex_tree_node_explicit_storage(Siblings* sib = nullptr, const Filtration_value& filtration = Filtration_value())
      : SimplexTree::Filtration_simplex_base(filtration), children_(sib)
  {}

  //will fail to compile when called with SimplexTree::Options::store_key is false.
  Simplex_tree_node_explicit_storage(Siblings* sib, Filtration_value filtration, Simplex_key key)
      : SimplexTree::Filtration_simplex_base(filtration), SimplexTree::Key_simplex_base(key), children_(sib)
  {}

  /*
   * Assign children to the node
   */
  void assign_children(Siblings * children) {
    children_ = children;
  }

  /* Careful -> children_ can be NULL*/
  Siblings * children() {
    return children_;
  }

  Simplex_data& data() { return boost::empty_value<Simplex_data>::get(); }

 private:
  Siblings * children_;
};

/** @}*/  // end addtogroup simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_
