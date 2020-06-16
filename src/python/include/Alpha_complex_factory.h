/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_ALPHA_COMPLEX_FACTORY_H_
#define INCLUDE_ALPHA_COMPLEX_FACTORY_H_

#include <gudhi/Simplex_tree.h>
#include <gudhi/Alpha_complex.h>
#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Alpha_complex_options.h>
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

template <typename CgalPointType>
std::vector<double> pt_cgal_to_cython(CgalPointType const& point) {
  std::vector<double> vd;
  for (auto coord = point.cartesian_begin(); coord != point.cartesian_end(); coord++)
    vd.push_back(CGAL::to_double(*coord));
  return vd;
}

template <typename CgalPointType>
static CgalPointType pt_cython_to_cgal(std::vector<double> const& vec) {
  return CgalPointType(vec.size(), vec.begin(), vec.end());
}

class Abstract_alpha_complex {
 public:
  virtual std::vector<double> get_point(int vh) = 0;
  virtual void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                                   bool default_filtration_value) = 0;
};

class Exact_Alphacomplex_dD : public Abstract_alpha_complex {
 private:
  using Exact_kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
  using Point_exact_kernel = typename Exact_kernel::Point_d;

 public:
  Exact_Alphacomplex_dD(const std::vector<std::vector<double>>& points, bool exact_version)
      : exact_version_(exact_version) {
    ac_exact_ptr_ = std::make_unique<Alpha_complex<Exact_kernel>>(
      boost::adaptors::transform(points, pt_cython_to_cgal<Point_exact_kernel>));
  }

  std::vector<double> get_point(int vh) {
    Point_exact_kernel const& point = ac_exact_ptr_->get_point(vh);
    return pt_cgal_to_cython(point);
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                           bool default_filtration_value){
    ac_exact_ptr_->create_complex(*simplex_tree, max_alpha_square, exact_version_, default_filtration_value);
  }

 private:
  bool exact_version_;
  std::unique_ptr<Alpha_complex<Exact_kernel>> ac_exact_ptr_;
};

class Inexact_Alphacomplex_dD : public Abstract_alpha_complex {
 private:
  using Inexact_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Point_inexact_kernel = typename Inexact_kernel::Point_d;

 public:
  Inexact_Alphacomplex_dD(const std::vector<std::vector<double>>& points, bool exact_version)
      : exact_version_(exact_version) {
    ac_inexact_ptr_ = std::make_unique<Alpha_complex<Inexact_kernel>>(
      boost::adaptors::transform(points, pt_cython_to_cgal<Point_inexact_kernel>));
  }

  std::vector<double> get_point(int vh) {
    Point_inexact_kernel const& point = ac_inexact_ptr_->get_point(vh);
    return pt_cgal_to_cython(point);
  }
  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                           bool default_filtration_value){
    ac_inexact_ptr_->create_complex(*simplex_tree, max_alpha_square, exact_version_, default_filtration_value);
  }

 private:
  bool exact_version_;
  std::unique_ptr<Alpha_complex<Inexact_kernel>> ac_inexact_ptr_;
};

template <complexity Complexity>
class Alphacomplex_3D : public Abstract_alpha_complex {
 private:
  using Point_3 = typename Alpha_complex_3d<Complexity, false, false>::Bare_point_3;

  static Point_3 pt_cython_to_cgal_3(std::vector<double> const& vec) {
    return Point_3(vec[0], vec[1], vec[2]);
  }

 public:
  Alphacomplex_3D(const std::vector<std::vector<double>>& points) {
    alpha3d_ptr_ = std::make_unique<Alpha_complex_3d<Complexity, false, false>>(
      boost::adaptors::transform(points, pt_cython_to_cgal_3));
  }

  std::vector<double> get_point(int vh) {
    Point_3 const& point = alpha3d_ptr_->get_point(vh);
    return pt_cgal_to_cython(point);
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                           bool default_filtration_value){
    alpha3d_ptr_->create_complex(*simplex_tree, max_alpha_square);
    if (default_filtration_value) {
      // TODO
    }
  }

 private:
  std::unique_ptr<Alpha_complex_3d<Complexity, false, false>> alpha3d_ptr_;
};


}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // INCLUDE_ALPHA_COMPLEX_FACTORY_H_
