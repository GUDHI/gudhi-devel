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
auto pairify(void* p, ssize_t h, ssize_t w) {
  return [=](ssize_t i){
    char* birth = (char*)p + i * h;
    char* death = birth + w;
    return std::make_pair(*(double*)birth, *(double*)death);
  };
}

double wasserstein_distance(
    Dgm d1, Dgm d2,
    double wasserstein_power, double internal_p,
    double delta)
{
  py::buffer_info buf1 = d1.request();
  py::buffer_info buf2 = d2.request();

  py::gil_scoped_release release;

  // shape (n,2) or (0) for empty
  if((buf1.ndim!=2 || buf1.shape[1]!=2) && (buf1.ndim!=1 || buf1.shape[0]!=0))
    throw std::runtime_error("Diagram 1 must be an array of size n x 2");
  if((buf2.ndim!=2 || buf2.shape[1]!=2) && (buf2.ndim!=1 || buf2.shape[0]!=0))
    throw std::runtime_error("Diagram 2 must be an array of size n x 2");
  ssize_t stride11 = buf1.ndim == 2 ? buf1.strides[1] : 0;
  ssize_t stride21 = buf2.ndim == 2 ? buf2.strides[1] : 0;
  auto cnt1 = boost::counting_range<ssize_t>(0, buf1.shape[0]);
  auto diag1 = boost::adaptors::transform(cnt1, pairify(buf1.ptr, buf1.strides[0], stride11));
  auto cnt2 = boost::counting_range<ssize_t>(0, buf2.shape[0]);
  auto diag2 = boost::adaptors::transform(cnt2, pairify(buf2.ptr, buf2.strides[0], stride21));

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
