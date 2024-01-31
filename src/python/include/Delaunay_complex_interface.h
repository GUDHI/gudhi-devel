/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_DELAUNAY_COMPLEX_INTERFACE_H_
#define INCLUDE_DELAUNAY_COMPLEX_INTERFACE_H_

#include "Delaunay_complex_factory.h"
#include <gudhi/Alpha_complex_options.h>

#include "Simplex_tree_interface.h"

#include <iostream>
#include <vector>
#include <string>
#include <memory>  // for std::unique_ptr

namespace Gudhi {

namespace delaunay_complex {

class Delaunay_complex_interface {
 public:
  Delaunay_complex_interface(const std::vector<std::vector<double>>& points,
                          const std::vector<double>& weights,
                          bool fast_version, bool exact_version) {
    const bool weighted = (weights.size() > 0);
    if (fast_version) {
      if (weighted) {
        alpha_ptr_ = std::make_unique<Inexact_delaunay_complex_dD<true>>(points, weights);
      } else {
        alpha_ptr_ = std::make_unique<Inexact_delaunay_complex_dD<false>>(points);
      }
    } else {
      if (weighted) {
        alpha_ptr_ = std::make_unique<Exact_delaunay_complex_dD<true>>(points, weights, exact_version);
      } else {
        alpha_ptr_ = std::make_unique<Exact_delaunay_complex_dD<false>>(points, exact_version);
      }
    }
  }

  std::vector<double> get_point(int vh) {
    return alpha_ptr_->get_point(vh);
  }

  void create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square,
                           bool default_filtration_value) {
    // Nothing to be done in case of an empty point set
    if (alpha_ptr_->num_vertices() > 0)
      alpha_ptr_->create_simplex_tree(simplex_tree, max_alpha_square, default_filtration_value);
  }

  static void set_float_relative_precision(double precision) {
    // cf. Exact_delaunay_complex_dD kernel type in Alpha_complex_factory.h
    CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::FT::set_relative_precision_of_to_double(precision);
  }

  static double get_float_relative_precision() {
    // cf. Exact_delaunay_complex_dD kernel type in Alpha_complex_factory.h
    return CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::FT::get_relative_precision_of_to_double();
  }

 private:
  std::unique_ptr<Abstract_delaunay_complex> alpha_ptr_;
};

}  // namespace delaunay_complex

}  // namespace Gudhi

#endif  // INCLUDE_DELAUNAY_COMPLEX_INTERFACE_H_
