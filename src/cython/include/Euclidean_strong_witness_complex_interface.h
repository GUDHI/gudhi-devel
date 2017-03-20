/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef EUCLIDEAN_STRONG_WITNESS_COMPLEX_INTERFACE_H
#define	EUCLIDEAN_STRONG_WITNESS_COMPLEX_INTERFACE_H

#include <gudhi/Simplex_tree.h>
#include <gudhi/Euclidean_strong_witness_complex.h>

#include "Simplex_tree_interface.h"

#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>  // std::pair
#include <iostream>
#include <cstddef>

namespace Gudhi {

namespace witness_complex {


class Euclidean_strong_witness_complex_interface {
  using Dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
  using Point_d = Dynamic_kernel::Point_d;

  typedef typename Simplex_tree<>::Simplex_key Simplex_key;

 public:
  Euclidean_strong_witness_complex_interface(const std::vector<std::vector<double>>& landmarks,
                                             const std::vector<std::vector<double>>& witnesses)
    : landmarks_(landmarks.begin(), landmarks.end()),
      witnesses_(witnesses.begin(), witnesses.end()),
      witness_complex_(landmarks_, witnesses_) {
  }

  void create_simplex_tree(Gudhi::Simplex_tree<>* simplex_tree, double max_alpha_square,
                           std::size_t limit_dimension) {
    witness_complex_.create_complex(*simplex_tree, max_alpha_square, limit_dimension);
    simplex_tree->initialize_filtration();
  }

  void create_simplex_tree(Gudhi::Simplex_tree<>* simplex_tree, double max_alpha_square) {
    witness_complex_.create_complex(*simplex_tree, max_alpha_square);
    simplex_tree->initialize_filtration();
  }

  std::vector<double> get_point(unsigned vh) {
    std::vector<double> vd;
    if (vh < landmarks_.size()) {
      Point_d ph = witness_complex_.get_point(vh);
      for (auto coord = ph.cartesian_begin(); coord < ph.cartesian_end(); coord++)
        vd.push_back(*coord);
    }
    return vd;
  }

 private:
  std::vector<Point_d> landmarks_;
  std::vector<Point_d> witnesses_;
  Euclidean_strong_witness_complex<Dynamic_kernel> witness_complex_;
};

}  // namespace witness_complex

} // namespace Gudhi

#endif  // EUCLIDEAN_STRONG_WITNESS_COMPLEX_INTERFACE_H

