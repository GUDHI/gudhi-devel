/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <wasserstein.h> // Hera

#include <utility>

namespace py = pybind11;
typedef py::array_t<double> Dgm;

// Get m[i,0] and m[i,1] as a pair
static auto pairify(void* p, ssize_t h, ssize_t w) {
  return [=](ssize_t i){
    char* birth = (char*)p + i * h;
    char* death = birth + w;
    return std::make_pair(*(double*)birth, *(double*)death);
  };
}

inline auto numpy_to_range_of_pairs(py::array_t<double> dgm) {
  py::buffer_info buf = dgm.request();
  // shape (n,2) or (0) for empty
  if((buf.ndim!=2 || buf.shape[1]!=2) && (buf.ndim!=1 || buf.shape[0]!=0))
    throw std::runtime_error("Diagram must be an array of size n x 2");
  // In the case of shape (0), avoid reading non-existing strides[1] even if we won't use it.
  ssize_t stride1 = buf.ndim == 2 ? buf.strides[1] : 0;
  auto cnt = boost::counting_range<ssize_t>(0, buf.shape[0]);
  return boost::adaptors::transform(cnt, pairify(buf.ptr, buf.strides[0], stride1));
  // Be careful that the returned range cannot contain references to dead temporaries.
}

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

PYBIND11_MODULE(hera, m) {
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
