/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_ALPHA_COMPLEX_INTERFACE_H_
#define INCLUDE_ALPHA_COMPLEX_INTERFACE_H_

#include "Alpha_complex_factory.h"
#include <gudhi/Alpha_complex_options.h>

#include "Simplex_tree_interface.h"

#include <iostream>
#include <vector>
#include <string>
#include <memory>  // for std::unique_ptr

namespace Gudhi {

namespace alpha_complex {

class Alpha_complex_interface {
 public:
  Alpha_complex_interface(const std::vector<std::vector<double>>& points, bool fast_version, bool exact_version) {
    if (points[0].size() == 3) {
      if (fast_version)
        alpha_ptr_ = std::make_unique<Alphacomplex_3D<Gudhi::alpha_complex::complexity::FAST>>(points);
      else if (exact_version)
        alpha_ptr_ = std::make_unique<Alphacomplex_3D<Gudhi::alpha_complex::complexity::EXACT>>(points);
      else
        alpha_ptr_ = std::make_unique<Alphacomplex_3D<Gudhi::alpha_complex::complexity::SAFE>>(points);
    } else {
      if (fast_version) {
        alpha_ptr_ = std::make_unique<Inexact_Alphacomplex_dD>(points, exact_version);
      } else {
        alpha_ptr_ = std::make_unique<Exact_Alphacomplex_dD>(points, exact_version);
      }
    }
  }

  std::vector<double> get_point(int vh) {
    return alpha_ptr_->get_point(vh);
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                           bool default_filtration_value) {
      alpha_ptr_->create_simplex_tree(simplex_tree, max_alpha_square, default_filtration_value);
  }

 private:
  std::unique_ptr<Abstract_alpha_complex> alpha_ptr_;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // INCLUDE_ALPHA_COMPLEX_INTERFACE_H_
