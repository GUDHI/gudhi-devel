#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <boost/range/iterator_range.hpp>

#include <wasserstein.h>

#include <array>

namespace py = pybind11;
typedef py::array_t<double, py::array::c_style | py::array::forcecast> Dgm;

double wasserstein_distance(
    Dgm d1, Dgm d2,
    double wasserstein_power, double internal_p,
    double delta)
{
  py::buffer_info buf1 = d1.request();
  py::buffer_info buf2 = d2.request();
  // shape (n,2) or (0) for empty
  if((buf1.ndim!=2 || buf1.shape[1]!=2) && (buf1.ndim!=1 || buf1.shape[0]!=0))
    throw std::runtime_error("Diagram 1 must be an array of size n x 2");
  if((buf2.ndim!=2 || buf2.shape[1]!=2) && (buf2.ndim!=1 || buf2.shape[0]!=0))
    throw std::runtime_error("Diagram 2 must be an array of size n x 2");
  typedef std::array<double, 2> Point;
  auto p1 = (Point*)buf1.ptr;
  auto p2 = (Point*)buf2.ptr;
  auto diag1 = boost::make_iterator_range(p1, p1+buf1.shape[0]);
  auto diag2 = boost::make_iterator_range(p2, p2+buf2.shape[0]);

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
          // Should we name those q, p and d instead?
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
