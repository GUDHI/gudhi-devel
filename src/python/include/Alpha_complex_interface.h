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

#include <gudhi/Simplex_tree.h>
#include <gudhi/Alpha_complex.h>
#include <CGAL/Epeck_d.h>
#include <CGAL/Epick_d.h>

#include <boost/range/adaptor/transformed.hpp>

#include "Simplex_tree_interface.h"

#include <iostream>
#include <vector>
#include <string>
#include <memory>  // for std::unique_ptr

namespace Gudhi {

namespace alpha_complex {

using Exact_kernel = CGAL::Epeck_d< CGAL::Dynamic_dimension_tag >;
using Inexact_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;

template <class Kernel>
class Alpha_complex_interface {
  using Point_d = typename Kernel::Point_d;

 public:
  Alpha_complex_interface(const std::vector<std::vector<double>>& points) {
    auto mkpt = [](std::vector<double> const& vec){
      return Point_d(vec.size(), vec.begin(), vec.end());
    };
    alpha_complex_ = std::make_unique<Alpha_complex<Kernel>>(boost::adaptors::transform(points, mkpt));
  }

  Alpha_complex_interface(const std::string& off_file_name, bool from_file = true) {
    alpha_complex_ = std::make_unique<Alpha_complex<Kernel>>(off_file_name);
  }

  std::vector<double> get_point(int vh) {
    std::vector<double> vd;
    Point_d const& ph = alpha_complex_->get_point(vh);
    for (auto coord = ph.cartesian_begin(); coord != ph.cartesian_end(); coord++)
      vd.push_back(CGAL::to_double(*coord));
    return vd;
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                           bool exact_version, bool default_filtration_value) {
    alpha_complex_->create_complex(*simplex_tree, max_alpha_square, exact_version, default_filtration_value);
  }

 private:
  std::unique_ptr<Alpha_complex<Kernel>> alpha_complex_;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // INCLUDE_ALPHA_COMPLEX_INTERFACE_H_
