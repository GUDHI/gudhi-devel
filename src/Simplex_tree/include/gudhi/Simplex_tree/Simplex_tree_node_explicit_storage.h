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

#ifndef SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_
#define SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_

#include <vector>

namespace Gudhi {

/* \addtogroup simplex_tree
 * Represents a node of a Simplex_tree.
 * @{
 */

/*
 * \brief Node of a simplex tree with filtration value
 * and simplex key.
 *
 * It stores explicitely its own filtration value and its own Simplex_key.
 */
template<class SimplexTree>
class Simplex_tree_node_explicit_storage {
 public:
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Filtration_value Filtration_value;
  typedef typename SimplexTree::Simplex_key Simplex_key;

  // Default constructor.
  Simplex_tree_node_explicit_storage()
      : children_(NULL),
        simplex_key_(-1),
        filtration_(0) {
  }

  Simplex_tree_node_explicit_storage(Siblings * sib,
                                     Filtration_value filtration)
      : children_(sib),
        simplex_key_(-1),
        filtration_(filtration) {
  }

  void assign_key(Simplex_key key) {
    simplex_key_ = key;
  }

  /*
   * Assign a children to the node
   */
  void assign_children(Siblings * children) {
    children_ = children;
  }
  /*
   *
   */
  void assign_filtration(double filtration_value) {
    filtration_ = filtration_value;
  }

  Filtration_value filtration() {
    return filtration_;
  }

  /* Careful -> children_ can be NULL*/
  Siblings * children() {
    return children_;
  }

  Simplex_key key() {
    return simplex_key_;
  }

 private:
  Siblings * children_;

  // Data attached to simplex, explicit storage
  Simplex_key simplex_key_;
  Filtration_value filtration_;   // value in the filtration
};

/* @} */  // end addtogroup simplex_tree
}  // namespace Gudhi

#endif  // SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_
