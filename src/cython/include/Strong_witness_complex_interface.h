/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
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

#ifndef INCLUDE_STRONG_WITNESS_COMPLEX_INTERFACE_H_
#define INCLUDE_STRONG_WITNESS_COMPLEX_INTERFACE_H_

#include <gudhi/Simplex_tree.h>
#include <gudhi/Strong_witness_complex.h>

#include "Simplex_tree_interface.h"

#include <vector>
#include <utility>  // std::pair
#include <iostream>
#include <cstddef>

namespace Gudhi {

namespace witness_complex {

class Strong_witness_complex_interface {
  using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
  using Nearest_landmark_table = std::vector<Nearest_landmark_range>;

 public:
  Strong_witness_complex_interface(const Nearest_landmark_table& nlt) {
    witness_complex_ = new Strong_witness_complex<Nearest_landmark_table>(nlt);
  }

  ~Strong_witness_complex_interface() {
    delete witness_complex_;
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double  max_alpha_square,
                           std::size_t limit_dimension) {
    witness_complex_->create_complex(*simplex_tree, max_alpha_square, limit_dimension);
    simplex_tree->initialize_filtration();
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree,
                           double  max_alpha_square) {
    witness_complex_->create_complex(*simplex_tree, max_alpha_square);
    simplex_tree->initialize_filtration();
  }

 private:
  Strong_witness_complex<Nearest_landmark_table>* witness_complex_;
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // INCLUDE_STRONG_WITNESS_COMPLEX_INTERFACE_H_
