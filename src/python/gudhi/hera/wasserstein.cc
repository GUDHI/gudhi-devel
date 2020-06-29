/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <wasserstein.h> // Hera

#include <pybind11_diagram_utils.h>

double wasserstein_distance(
    Dgm d1, Dgm d2,
    double wasserstein_power, double internal_p,
    double delta)
{
  // I *think* the call to request() has to be before releasing the GIL.
  auto diag1 = numpy_to_range_of_pairs(d1);
  auto diag2 = numpy_to_range_of_pairs(d2);

  py::gil_scoped_release release;

  hera::AuctionParams<double> params;
  params.wasserstein_power = wasserstein_power;
  // hera encodes infinity as -1...
  if(std::isinf(internal_p)) internal_p = hera::get_infinity<double>();
  params.internal_p = internal_p;
  params.delta = delta;
  // The extra parameters are purposedly not exposed for now.
  return hera::wasserstein_dist(diag1, diag2, params);
}

PYBIND11_MODULE(wasserstein, m) {
      m.def("wasserstein_distance", &wasserstein_distance,
          py::arg("X"), py::arg("Y"),
          py::arg("order") = 1,
          py::arg("internal_p") = std::numeric_limits<double>::infinity(),
          py::arg("delta") = .01,
          R"pbdoc(
        Compute the Wasserstein distance between two diagrams.
        Points at infinity are supported.

        Parameters:
            X (n x 2 numpy array): First diagram
            Y (n x 2 numpy array): Second diagram
            order (float): Wasserstein exponent W_q
            internal_p (float): Internal Minkowski norm L^p in R^2
            delta (float): Relative error 1+delta

        Returns:
            float: Approximate Wasserstein distance W_q(X,Y)
    )pbdoc");
}
