/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2025/01 Vincent Rouvreau: Use nanobind instead of PyBind11 for python bindings
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Bottleneck.h>
//#include <optional>
#include <pybind11_diagram_utils.h>
#include <nanobind/stl/optional.h>
//#include <pybind11/stl.h>

// Indices are added internally in bottleneck_distance, they are not needed in the input.
static auto make_point(double x, double y, nb::ssize_t) { return std::pair(x, y); };

// For compatibility with older versions, we want to support e=None.
double bottleneck(const Dgm& d1, const Dgm& d2, std::optional<double> epsilon)
{
  double e = epsilon.value_or((std::numeric_limits<double>::min)());
// I *think* the call to request() in numpy_to_range_of_pairs has to be before releasing the GIL.
  auto diag1 = numpy_to_range_of_pairs(d1, make_point);
  auto diag2 = numpy_to_range_of_pairs(d2, make_point);

  nb::gil_scoped_release release;

  return Gudhi::persistence_diagram::bottleneck_distance(diag1, diag2, e);
}

NB_MODULE(bottleneck, m) {
      m.attr("__license__") = "GPL v3";
      m.def("bottleneck_distance", &bottleneck,
          nb::arg("diagram_1"), nb::arg("diagram_2"),
          nb::arg("e") = nb::none(),
          R"pbdoc(
    Compute the Bottleneck distance between two diagrams.
    Points at infinity and on the diagonal are supported.

    :param diagram_1: The first diagram.
    :type diagram_1: numpy array of shape (m,2)
    :param diagram_2: The second diagram.
    :type diagram_2: numpy array of shape (n,2)
    :param e: If `e` is 0, this uses an expensive algorithm to compute the
        exact distance.
        If `e` is not 0, it asks for an additive `e`-approximation, and
        currently also allows a small multiplicative error (the last 2 or 3
        bits of the mantissa may be wrong). This version of the algorithm takes
        advantage of the limited precision of `double` and is usually a lot
        faster to compute, whatever the value of `e`.
        Thus, by default (`e=None`), `e` is the smallest positive double.
    :type e: float
    :rtype: float
    :returns: the bottleneck distance.
    )pbdoc");
}
