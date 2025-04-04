/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2025/01 Vincent Rouvreau: Use nanobind instead of PyBind11 for python bindings
 *      - 2025/04 Hannah Schreiber: Re-add possibility of native Python sequences as input
 *      - YYYY/MM Author: Description of the modification
 */

#include <limits>
#include <optional>
#include <utility>  //std::pair

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>

#include <gudhi/Bottleneck.h>
#include <python_interfaces/diagram_utils.h>

namespace nb = nanobind;

// Indices are added internally in bottleneck_distance, they are not needed in the input.
static auto _make_point(double x, double y, std::size_t) { return std::pair(x, y); };

// For compatibility with older versions, we want to support e=None.
template <class Dgm>
double _bottleneck(const Dgm& d1, const Dgm& d2, std::optional<double> epsilon)
{
  double e = epsilon.value_or((std::numeric_limits<double>::min)());
  // I *think* the call to request() in array_to_range_of_pairs has to be before releasing the GIL.
  auto diag1 = array_to_range_of_pairs(d1, _make_point);
  auto diag2 = array_to_range_of_pairs(d2, _make_point);

  nb::gil_scoped_release release;

  return Gudhi::persistence_diagram::bottleneck_distance(diag1, diag2, e);
}

NB_MODULE(_bottleneck_ext, m)
{
  m.attr("__license__") = "GPL v3";
  m.def("_bottleneck_distance_tensor",
        &_bottleneck<Tensor_dgm>,
        nb::arg("diagram_1"),
        nb::arg("diagram_2"),
        nb::arg("e") = nb::none())
      .def("_bottleneck_distance_list",
           &_bottleneck<List_dgm>,
           nb::arg("diagram_1"),
           nb::arg("diagram_2"),
           nb::arg("e") = nb::none())
      .def("_bottleneck_distance_sequence",
           &_bottleneck<Sequence_dgm>,
           nb::arg("diagram_1"),
           nb::arg("diagram_2"),
           nb::arg("e") = nb::none());
}
