/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/03 Vincent Rouvreau: Use nanobind instead of Cython for python bindings.
 *      - YYYY/MM Author: Description of the modification
 */


#include <gudhi/Simplex_tree.h>
#include <gudhi/Euclidean_witness_complex.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/ndarray.h>

#include "Simplex_tree_interface_cython.h"

#include <CGAL/Epick_d.h>

#include <vector>
#include <iostream>
#include <cstddef>
#include <limits>

namespace Gudhi {

namespace witness_complex {


class Euclidean_witness_complex_interface {
  using Dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
  using Point_d = Dynamic_kernel::Point_d;

 public:
  Euclidean_witness_complex_interface(const std::vector<std::vector<double>>& landmarks,
                                      const std::vector<std::vector<double>>& witnesses) {
    landmarks_.reserve(landmarks.size());
    for (auto& landmark : landmarks)
      landmarks_.emplace_back(landmark.begin(), landmark.end());
    witness_complex_ = std::make_unique<Euclidean_witness_complex<Dynamic_kernel>>(landmarks_, witnesses);
  }

  void create_simplex_tree(Simplex_tree_interface& simplex_tree,
                           double max_alpha_square,
                           std::size_t limit_dimension = std::numeric_limits<std::size_t>::max()) {
    witness_complex_->create_complex(simplex_tree, max_alpha_square, limit_dimension);
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
  std::unique_ptr<Euclidean_witness_complex<Dynamic_kernel>> witness_complex_;
};

}  // namespace witness_complex

}  // namespace Gudhi


namespace nb = nanobind;

using Kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
using Point_d = Kernel::Point_d;

namespace egwc = Gudhi::witness_complex;
using egwci = egwc::Euclidean_witness_complex_interface;


NB_MODULE(_euclidean_witness_complex_ext, m) {
    m.attr("__license__") = "GPL v3";

    nb::class_<egwci>(m, "Euclidean_witness_complex_interface")
        .def(nb::init<const std::vector<std::vector<double>>&,
                      const std::vector<std::vector<double>>&>(), "Constructor")
        .def("create_simplex_tree", &egwci::create_simplex_tree, "")
        .def("get_point", &egwci::get_point, "");
}
