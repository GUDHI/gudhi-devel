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
double _bottleneck(const Tensor_dgm& d1, const Tensor_dgm& d2, std::optional<double> epsilon)
{
  double e = epsilon.value_or((std::numeric_limits<double>::min)());
  auto d1_view = d1.view(); // views have to live until the end of the method, as not copied in
  auto d2_view = d2.view(); // the iterator of array_to_range_of_pairs
  // I *think* the call to request() in array_to_range_of_pairs has to be before releasing the GIL.
  auto diag1 = array_to_range_of_pairs(d1_view, _make_point);
  auto diag2 = array_to_range_of_pairs(d2_view, _make_point);

  nb::gil_scoped_release release;

  return Gudhi::persistence_diagram::bottleneck_distance(diag1, diag2, e);
}

NB_MODULE(_bottleneck_ext, m)
{
  m.attr("__license__") = "GPL v3";
  m.def("_bottleneck_distance",
        &_bottleneck,
        nb::arg("diagram_1"),
        nb::arg("diagram_2"),
        nb::arg("e") = nb::none());
}
