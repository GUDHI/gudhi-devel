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

// template Functor that transforms a CGAL point to a vector of double as expected by cython
template<typename CgalPointType, bool Weighted>
struct Point_cgal_to_cython;

// Specialized Unweighted Functor
template<typename CgalPointType>
struct Point_cgal_to_cython<CgalPointType, false> {
  std::vector<double> operator()(CgalPointType const& point) const 
  {
    std::vector<double> vd;
    vd.reserve(point.dimension());
    for (auto coord = point.cartesian_begin(); coord != point.cartesian_end(); coord++)
      vd.push_back(CGAL::to_double(*coord));
    return vd;
  }
};

// Specialized Weighted Functor
template<typename CgalPointType>
struct Point_cgal_to_cython<CgalPointType, true> {
  std::vector<double> operator()(CgalPointType const& weighted_point) const 
  {
    const auto& point = weighted_point.point();
    return Point_cgal_to_cython<decltype(point), false>()(point);
  }
};

// Function that transforms a cython point (aka. a vector of double) to a CGAL point
template <typename CgalPointType>
static CgalPointType pt_cython_to_cgal(std::vector<double> const& vec) {
  return CgalPointType(vec.size(), vec.begin(), vec.end());
}

class Abstract_alpha_complex {
 public:
  virtual std::vector<double> get_point(int vh) = 0;

  virtual bool create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square,
                                   bool default_filtration_value) = 0;
  
  virtual std::size_t num_vertices() const = 0;

  virtual ~Abstract_alpha_complex() = default;
};

template <bool Weighted = false>
class Exact_alpha_complex_dD final : public Abstract_alpha_complex {
 private:
  using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
  using Bare_point = typename Kernel::Point_d;
  using Point = std::conditional_t<Weighted, typename Kernel::Weighted_point_d,
                                             typename Kernel::Point_d>;

 public:
  Exact_alpha_complex_dD(const std::vector<std::vector<double>>& points, bool exact_version)
    : exact_version_(exact_version),
      alpha_complex_(boost::adaptors::transform(points, pt_cython_to_cgal<Bare_point>)) {
  }

  Exact_alpha_complex_dD(const std::vector<std::vector<double>>& points,
                           const std::vector<double>& weights, bool exact_version)
    : exact_version_(exact_version),
      alpha_complex_(boost::adaptors::transform(points, pt_cython_to_cgal<Bare_point>), weights) {
  }

  virtual std::vector<double> get_point(int vh) override {
    // Can be a Weighted or a Bare point in function of Weighted
    return Point_cgal_to_cython<Point, Weighted>()(alpha_complex_.get_point(vh));
  }

  virtual bool create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square,
                                   bool default_filtration_value) override {
    return alpha_complex_.create_complex(*simplex_tree, max_alpha_square, exact_version_, default_filtration_value);
  }

  virtual std::size_t num_vertices() const override {
    return alpha_complex_.num_vertices();
  }

 private:
  bool exact_version_;
  Alpha_complex<Kernel, Weighted> alpha_complex_;
};

template <bool Weighted = false>
class Inexact_alpha_complex_dD final : public Abstract_alpha_complex {
 private:
  using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Bare_point = typename Kernel::Point_d;
  using Point = std::conditional_t<Weighted, typename Kernel::Weighted_point_d,
                                             typename Kernel::Point_d>;

 public:
  Inexact_alpha_complex_dD(const std::vector<std::vector<double>>& points)
    : alpha_complex_(boost::adaptors::transform(points, pt_cython_to_cgal<Bare_point>)) {
  }

  Inexact_alpha_complex_dD(const std::vector<std::vector<double>>& points, const std::vector<double>& weights)
    : alpha_complex_(boost::adaptors::transform(points, pt_cython_to_cgal<Bare_point>), weights) {
  }

  virtual std::vector<double> get_point(int vh) override {
    // Can be a Weighted or a Bare point in function of Weighted
    return Point_cgal_to_cython<Point, Weighted>()(alpha_complex_.get_point(vh));
  }
  virtual bool create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square,
                                   bool default_filtration_value) override {
    return alpha_complex_.create_complex(*simplex_tree, max_alpha_square, false, default_filtration_value);
  }

  virtual std::size_t num_vertices() const override {
    return alpha_complex_.num_vertices();
  }

 private:
  Alpha_complex<Kernel, Weighted> alpha_complex_;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // INCLUDE_ALPHA_COMPLEX_FACTORY_H_
