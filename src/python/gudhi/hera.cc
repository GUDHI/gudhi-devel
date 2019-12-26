#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <boost/range/iterator_range.hpp>

#include <wasserstein.h>

#include <array>

namespace py = pybind11;
typedef py::array_t<double, py::array::c_style | py::array::forcecast> Dgm;

namespace hera {
template <> struct DiagramTraits<Dgm>{
  //using Container = void;
  using PointType = std::array<double,2>;
  using RealType  = double;

  static RealType get_x(const PointType& p) { return std::get<0>(p); }
  static RealType get_y(const PointType& p) { return std::get<1>(p); }
};
}

double wasserstein_distance(
    Dgm d1,
    Dgm d2)
{
  py::buffer_info buf1 = d1.request();
  py::buffer_info buf2 = d2.request();
  if(buf1.ndim!=2 || buf1.shape[1]!=2)
    throw std::runtime_error("Diagram 1 must be an array of size n x 2");
  if(buf2.ndim!=2 || buf2.shape[1]!=2)
    throw std::runtime_error("Diagram 1 must be an array of size n x 2");
  typedef hera::DiagramTraits<Dgm>::PointType Point;
  auto p1 = (Point*)buf1.ptr;
  auto p2 = (Point*)buf2.ptr;
  auto diag1 = boost::make_iterator_range(p1, p1+buf1.shape[0]);
  auto diag2 = boost::make_iterator_range(p2, p2+buf2.shape[0]);

  hera::AuctionParams<double> params;
  return hera::wasserstein_dist(diag1, diag2, params);
}

PYBIND11_MODULE(hera, m) {
      m.def("wasserstein_distance", &wasserstein_distance, R"pbdoc(
        Compute the Wasserstein distance between two diagrams
    )pbdoc");
}
