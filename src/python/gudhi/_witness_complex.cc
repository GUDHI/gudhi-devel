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

#include <vector>
#include <utility>  // std::pair
#include <cstddef>
#include <limits>  // for std::numeric_limits

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/ndarray.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <python_interfaces/Simplex_tree_interface.h>

using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
using Nearest_landmark_table = std::vector<Nearest_landmark_range>;

namespace nb = nanobind;
namespace gwc = Gudhi::witness_complex;
using gwci = gwc::Witness_complex<Nearest_landmark_table>;

std::size_t max_size_t = std::numeric_limits<std::size_t>::max();

NB_MODULE(_witness_complex_ext, m)
{
  m.attr("__license__") = "GPL v3";

  nb::class_<gwci>(m, "Witness_complex_interface")
      .def(nb::init<const Nearest_landmark_table&>(), "Constructor")
      .def(nb::init<>(), "Constructor")
      .def("create_simplex_tree",
           &gwci::create_complex<Gudhi::Simplex_tree_interface>,
           nb::arg("complex"),
           nb::arg("max_alpha_square"),
           nb::arg("limit_dimension") = max_size_t,
           "");
}
