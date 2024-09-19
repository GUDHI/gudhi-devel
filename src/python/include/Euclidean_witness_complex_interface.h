/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_EUCLIDEAN_WITNESS_COMPLEX_INTERFACE_H_
#define INCLUDE_EUCLIDEAN_WITNESS_COMPLEX_INTERFACE_H_

#include <gudhi/Simplex_tree.h>
#include <gudhi/Euclidean_witness_complex.h>

#include "Simplex_tree_interface.h"

#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>  // std::pair
#include <iostream>
#include <cstddef>

namespace Gudhi {

namespace witness_complex {


class Euclidean_witness_complex_interface {
  using Dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
  using Point_d = Dynamic_kernel::Point_d;

  typedef typename Simplex_tree<>::Simplex_key Simplex_key;

 public:
  Euclidean_witness_complex_interface(const std::vector<std::vector<double>>& landmarks,
                                      const std::vector<std::vector<double>>& witnesses) {
    landmarks_.reserve(landmarks.size());
    for (auto& landmark : landmarks)
      landmarks_.emplace_back(landmark.begin(), landmark.end());
    witness_complex_ = new Euclidean_witness_complex<Dynamic_kernel>(landmarks_, witnesses);
  }

  ~Euclidean_witness_complex_interface() {
    delete witness_complex_;
  }

  void create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square, std::size_t limit_dimension) {
    witness_complex_->create_complex(*simplex_tree, max_alpha_square, limit_dimension);
  }

  void create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square) {
    witness_complex_->create_complex(*simplex_tree, max_alpha_square);
  }

  std::vector<double> get_point(unsigned vh) {
    std::vector<double> vd;
    if (vh < landmarks_.size()) {
      Point_d ph = witness_complex_->get_point(vh);
      for (auto coord = ph.cartesian_begin(); coord < ph.cartesian_end(); coord++)
        vd.push_back(*coord);
    }
    return vd;
  }

 private:
  std::vector<Point_d> landmarks_;
  Euclidean_witness_complex<Dynamic_kernel>* witness_complex_;
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // INCLUDE_EUCLIDEAN_WITNESS_COMPLEX_INTERFACE_H_

