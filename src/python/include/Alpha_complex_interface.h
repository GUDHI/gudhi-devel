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

class Alpha_complex_interface {
 private:
  using Exact_kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
  using Inexact_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Point_exact_kernel = typename Exact_kernel::Point_d;
  using Point_inexact_kernel = typename Inexact_kernel::Point_d;

  template <typename CgalPointType>
  std::vector<double> pt_cgal_to_cython(CgalPointType& ph) {
    std::vector<double> vd;
    for (auto coord = ph.cartesian_begin(); coord != ph.cartesian_end(); coord++)
      vd.push_back(CGAL::to_double(*coord));
    return vd;
  }

  template <typename CgalPointType>
  static CgalPointType pt_cython_to_cgal(std::vector<double> const& vec) {
    return CgalPointType(vec.size(), vec.begin(), vec.end());
  }

 public:
  Alpha_complex_interface(const std::vector<std::vector<double>>& points, bool fast_version)
      : fast_version_(fast_version) {
    if (fast_version_) {
      ac_inexact_ptr_ = std::make_unique<Alpha_complex<Inexact_kernel>>(
          boost::adaptors::transform(points, pt_cython_to_cgal<Point_inexact_kernel>));
    } else {
      ac_exact_ptr_ = std::make_unique<Alpha_complex<Exact_kernel>>(
          boost::adaptors::transform(points, pt_cython_to_cgal<Point_exact_kernel>));
    }
  }

  Alpha_complex_interface(const std::string& off_file_name, bool fast_version, bool from_file = true)
      : fast_version_(fast_version) {
    if (fast_version_)
      ac_inexact_ptr_ = std::make_unique<Alpha_complex<Inexact_kernel>>(off_file_name);
    else
      ac_exact_ptr_ = std::make_unique<Alpha_complex<Exact_kernel>>(off_file_name);
  }

  std::vector<double> get_point(int vh) {
    if (fast_version_) {
      Point_inexact_kernel const& ph = ac_inexact_ptr_->get_point(vh);
      return pt_cgal_to_cython(ph);
    } else {
      Point_exact_kernel const& ph = ac_exact_ptr_->get_point(vh);
      return pt_cgal_to_cython(ph);
    }
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square, bool exact_version,
                           bool default_filtration_value) {
    if (fast_version_)
      ac_inexact_ptr_->create_complex(*simplex_tree, max_alpha_square, exact_version, default_filtration_value);
    else
      ac_exact_ptr_->create_complex(*simplex_tree, max_alpha_square, exact_version, default_filtration_value);
  }

 private:
  bool fast_version_;
  std::unique_ptr<Alpha_complex<Exact_kernel>> ac_exact_ptr_;
  std::unique_ptr<Alpha_complex<Inexact_kernel>> ac_inexact_ptr_;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // INCLUDE_ALPHA_COMPLEX_INTERFACE_H_
