/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/03 Vincent Rouvreau: Use nanobind instead of Cython for python bindings.
 *      - 2025/04 Hannah Schreiber: Re-add possibility of tensors (numpy, torch etc.) as input.
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

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <python_interfaces/Simplex_tree_interface.h>
#include <python_interfaces/points_utils.h>

namespace Gudhi {
namespace witness_complex {

class Witness_complex_interface
{
 public:
  using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
  using Nearest_landmark_table = std::vector<Nearest_landmark_range>;

 public:
  Witness_complex_interface() {}

  Witness_complex_interface(const Nearest_landmark_sequence& nlt) { witness_complex_.emplace(nlt); }

  // TODO: remove one copy by directly constructing std::vector<Point_d> instead of std::vector<std::vector<double>>
  Witness_complex_interface(const Nearest_landmark_tensor& nlt)
      : Witness_complex_interface(_get_sequence_from_tensor(nlt))
  {}

  void create_simplex_tree(Simplex_tree_interface* simplex_tree,
                           double max_alpha_square,
                           std::size_t limit_dimension = std::numeric_limits<std::size_t>::max())
  {
    if (witness_complex_) witness_complex_->create_complex(*simplex_tree, max_alpha_square, limit_dimension);
  }

 private:
  std::optional<Witness_complex<Nearest_landmark_table>> witness_complex_;
};

}  // namespace witness_complex
}  // namespace Gudhi

namespace nb = nanobind;
namespace gwc = Gudhi::witness_complex;
using gwci = gwc::Witness_complex_interface;

NB_MODULE(_witness_complex_ext, m)
{
  m.attr("__license__") = "MIT";

  nb::class_<gwci>(m, "Witness_complex_interface")
      .def(nb::init<const Nearest_landmark_sequence&>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Nearest_landmark_tensor&>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<>(), nb::call_guard<nb::gil_scoped_release>())
      .def("create_simplex_tree",
           &gwci::create_simplex_tree,
           nb::arg("simplex_tree"),
           nb::arg("max_alpha_square"),
           nb::arg("limit_dimension") = std::numeric_limits<std::size_t>::max(),
           nb::call_guard<nb::gil_scoped_release>());
}
