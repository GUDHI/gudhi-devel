/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Bottleneck.h>

#include <pybind11_diagram_utils.h>

double bottleneck(Dgm d1, Dgm d2, double epsilon)
{
  // I *think* the call to request() has to be before releasing the GIL.
  auto diag1 = numpy_to_range_of_pairs(d1);
  auto diag2 = numpy_to_range_of_pairs(d2);

  py::gil_scoped_release release;

  return Gudhi::persistence_diagram::bottleneck_distance(diag1, diag2, epsilon);
}

PYBIND11_MODULE(bottleneck, m) {
      m.attr("__license__") = "GPL v3";
      m.def("bottleneck_distance", &bottleneck,
          py::arg("diagram_1"), py::arg("diagram_2"),
          py::arg("e") = (std::numeric_limits<double>::min)(),
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
        Thus, by default, `e` is the smallest positive double.
    :type e: float
    :rtype: float
    :returns: the bottleneck distance.
    )pbdoc");
}
