/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

// https://github.com/grey-narn/hera/issues/3
// ssize_t is a non-standard type (well, posix)
// BaseTsd.h provides SSIZE_T on windows, this one should be the same there.
#ifdef _MSC_VER
#include <cstddef>
typedef std::ptrdiff_t ssize_t;
#endif

#include <bottleneck.h> // Hera

#include <pybind11_diagram_utils.h>

double bottleneck_distance(Dgm d1, Dgm d2, double delta)
{
  // I *think* the call to request() has to be before releasing the GIL.
  auto diag1 = numpy_to_range_of_pairs(d1);
  auto diag2 = numpy_to_range_of_pairs(d2);

  py::gil_scoped_release release;

  if (delta == 0)
    return hera::bottleneckDistExact(diag1, diag2);
  else
    return hera::bottleneckDistApprox(diag1, diag2, delta);
}

PYBIND11_MODULE(bottleneck, m) {
      m.def("bottleneck_distance", &bottleneck_distance,
          py::arg("X"), py::arg("Y"),
          py::arg("delta") = .01,
          R"pbdoc(
        Compute the Bottleneck distance between two diagrams.
        Points at infinity are supported.

        Parameters:
            X (n x 2 numpy array): First diagram
            Y (n x 2 numpy array): Second diagram
            delta (float): Relative error 1+delta

        Returns:
            float: (approximate) bottleneck distance d_B(X,Y)
    )pbdoc");
}
