/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2024/03 Vincent Rouvreau: Renamed Alpha_complex_interface as Delaunay_complex_interface for DelaunayCechComplex.
 *      - 2024/10 Vincent Rouvreau: Add square root filtration values interface
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_DELAUNAY_COMPLEX_INTERFACE_H_
#define INCLUDE_DELAUNAY_COMPLEX_INTERFACE_H_

#include "Delaunay_complex_factory.h"
#include <gudhi/Alpha_complex_options.h>
#include <gudhi/MEB_filtration.h>

#include "Simplex_tree_interface.h"

#include <CGAL/Epeck_d.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <memory>  // for std::unique_ptr
#include <cstddef>  // for std::size_t

namespace Gudhi {

namespace delaunay_complex {

class Delaunay_complex_interface {
 public:
  Delaunay_complex_interface(const std::vector<std::vector<double>>& points,
                             const std::vector<double>& weights,
                             bool fast_version, bool exact_version) {
    const bool weighted = (weights.size() > 0);
    // Specific cases for dimensions 2 and 3
    const std::size_t dimension = ((points.size() > 0) ? points[0].size() : 0);
    if (fast_version) {
      if (weighted) {
        if (dimension == 2) {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<2>>, true>>(points, weights, exact_version);
        } else if (dimension == 3) {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<3>>, true>>(points, weights, exact_version);
        } else {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>, true>>(points, weights, exact_version);
        }
      } else {
        if (dimension == 2) {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<2>>, false>>(points, exact_version);
        } else if (dimension == 3) {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<3>>, false>>(points, exact_version);
        } else {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>, false>>(points, exact_version);
        }
      }
    } else {
      if (weighted) {
        if (dimension == 2) {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<2>>, true>>(points, weights, exact_version);
        } else if (dimension == 3) {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<3>>, true>>(points, weights, exact_version);
        } else {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>, true>>(points, weights, exact_version);
        }
      } else {
        if (dimension == 2) {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<2>>, false>>(points, exact_version);
        } else if (dimension == 3) {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<3>>, false>>(points, exact_version);
        } else {
          delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>, false>>(points, exact_version);
        }
      }
    }
  }

  std::vector<double> get_point(int vh) {
    return delaunay_ptr_->get_point(vh);
  }

  void create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square,
                           Delaunay_filtration filtration, bool output_squared_values) {
    // Nothing to be done in case of an empty point set
    if (delaunay_ptr_->num_vertices() > 0)
      delaunay_ptr_->create_simplex_tree(simplex_tree, max_alpha_square, filtration, output_squared_values);
  }

  static void set_float_relative_precision(double precision) {
    // cf. CGAL::Epeck_d kernel type in Delaunay_complex_interface
    CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::FT::set_relative_precision_of_to_double(precision);
  }

  static double get_float_relative_precision() {
    // cf. CGAL::Epeck_d kernel type in Delaunay_complex_interface
    return CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::FT::get_relative_precision_of_to_double();
  }

 private:
  std::unique_ptr<Abstract_delaunay_complex> delaunay_ptr_;
};

}  // namespace delaunay_complex

}  // namespace Gudhi

#endif  // INCLUDE_DELAUNAY_COMPLEX_INTERFACE_H_
