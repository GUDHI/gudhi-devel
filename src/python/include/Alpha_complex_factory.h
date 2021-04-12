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
  vd.reserve(point.dimension());
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

  virtual bool create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                                   bool default_filtration_value) = 0;

  virtual ~Abstract_alpha_complex() = default;
};

class Exact_Alphacomplex_dD final : public Abstract_alpha_complex {
 private:
  using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
  using Point = typename Kernel::Point_d;

 public:
  Exact_Alphacomplex_dD(const std::vector<std::vector<double>>& points, bool exact_version)
    : exact_version_(exact_version),
      alpha_complex_(boost::adaptors::transform(points, pt_cython_to_cgal<Point>)) {
  }

  virtual std::vector<double> get_point(int vh) override {
    Point const& point = alpha_complex_.get_point(vh);
    return pt_cgal_to_cython(point);
  }

  virtual bool create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                                   bool default_filtration_value) override {
    return alpha_complex_.create_complex(*simplex_tree, max_alpha_square, exact_version_, default_filtration_value);
  }

 private:
  bool exact_version_;
  Alpha_complex<Kernel> alpha_complex_;
};

class Inexact_Alphacomplex_dD final : public Abstract_alpha_complex {
 private:
  using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Point = typename Kernel::Point_d;

 public:
  Inexact_Alphacomplex_dD(const std::vector<std::vector<double>>& points, bool exact_version)
    : exact_version_(exact_version),
      alpha_complex_(boost::adaptors::transform(points, pt_cython_to_cgal<Point>)) {
  }

  virtual std::vector<double> get_point(int vh) override {
    Point const& point = alpha_complex_.get_point(vh);
    return pt_cgal_to_cython(point);
  }
  virtual bool create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                                   bool default_filtration_value) override {
    return alpha_complex_.create_complex(*simplex_tree, max_alpha_square, exact_version_, default_filtration_value);
  }

 private:
  bool exact_version_;
  Alpha_complex<Kernel> alpha_complex_;
};

template <complexity Complexity>
class Alphacomplex_3D final : public Abstract_alpha_complex {
 private:
  using Point = typename Alpha_complex_3d<Complexity, false, false>::Bare_point_3;

  static Point pt_cython_to_cgal_3(std::vector<double> const& vec) {
    return Point(vec[0], vec[1], vec[2]);
  }

 public:
  Alphacomplex_3D(const std::vector<std::vector<double>>& points)
    : alpha_complex_(boost::adaptors::transform(points, pt_cython_to_cgal_3)) {
  }

  virtual std::vector<double> get_point(int vh) override {
    Point const& point = alpha_complex_.get_point(vh);
    return pt_cgal_to_cython(point);
  }

  virtual bool create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square,
                           bool default_filtration_value) override {
    return alpha_complex_.create_complex(*simplex_tree, max_alpha_square);
  }

 private:
  Alpha_complex_3d<Complexity, false, false> alpha_complex_;
};


}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // INCLUDE_ALPHA_COMPLEX_FACTORY_H_
