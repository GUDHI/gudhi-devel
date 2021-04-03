/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_ALPHA_COMPLEX_INTERFACE_3D_H_
#define INCLUDE_ALPHA_COMPLEX_INTERFACE_3D_H_

#include "Alpha_complex_factory.h"
#include <gudhi/Alpha_complex_options.h>

#include "Simplex_tree_interface.h"

#include <iostream>
#include <vector>
#include <string>
#include <memory>  // for std::unique_ptr

namespace Gudhi {

namespace alpha_complex {

class Alpha_complex_interface_3d {
 public:
  Alpha_complex_interface_3d(const std::vector<std::vector<double>>& points,
                          const std::vector<double>& weights,
                          bool fast_version, bool exact_version)
  : empty_point_set_(points.size() == 0) {
    const bool weighted = (weights.size() > 0);
    if (fast_version)
      if (weighted)
        alpha_ptr_ = std::make_unique<Alpha_complex_3D<Gudhi::alpha_complex::complexity::FAST, true>>(points, weights);
      else
        alpha_ptr_ = std::make_unique<Alpha_complex_3D<Gudhi::alpha_complex::complexity::FAST>>(points);
    else if (exact_version)
      if (weighted)
        alpha_ptr_ = std::make_unique<Alpha_complex_3D<Gudhi::alpha_complex::complexity::EXACT, true>>(points, weights);
      else
        alpha_ptr_ = std::make_unique<Alpha_complex_3D<Gudhi::alpha_complex::complexity::EXACT>>(points);
    else
      if (weighted)
        alpha_ptr_ = std::make_unique<Alpha_complex_3D<Gudhi::alpha_complex::complexity::SAFE, true>>(points, weights);
      else
        alpha_ptr_ = std::make_unique<Alpha_complex_3D<Gudhi::alpha_complex::complexity::SAFE>>(points);
  }

  std::vector<double> get_point(int vh) {
    return alpha_ptr_->get_point(vh);
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square) {
    // Nothing to be done in case of an empty point set
    if (!empty_point_set_)
      alpha_ptr_->create_simplex_tree(simplex_tree, max_alpha_square, false);
  }

 private:
  std::unique_ptr<Abstract_alpha_complex> alpha_ptr_;
  bool empty_point_set_;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // INCLUDE_ALPHA_COMPLEX_INTERFACE_3D_H_
