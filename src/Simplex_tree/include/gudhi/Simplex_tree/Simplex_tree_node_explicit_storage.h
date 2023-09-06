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

#include <vector>

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
template<class SimplexTree>
struct GUDHI_EMPTY_BASE_CLASS_OPTIMIZATION Simplex_tree_node_explicit_storage : SimplexTree::Filtration_simplex_base,
                                                                          SimplexTree::Key_simplex_base,
                                                                          SimplexTree::Hooks_simplex_base {
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Filtration_value Filtration_value;
  typedef typename SimplexTree::Simplex_key Simplex_key;

  Simplex_tree_node_explicit_storage(Siblings * sib = nullptr,
                                     const Filtration_value& filtration = {}) 
      : children_(sib) {
    this->assign_filtration(filtration);
  }

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

 private:
  Siblings * children_;
};

/** @}*/  // end addtogroup simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_
