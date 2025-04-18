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

namespace py = pybind11;
typedef py::array_t<double> Dgm;

// build_point(double birth, double death, ssize_t index) -> Point
template <class BuildPoint>
inline auto numpy_to_range_of_pairs(py::array_t<double> dgm, BuildPoint build_point)
{
  py::buffer_info buf = dgm.request();
  // shape (n,2) or (0) for empty
  if ((buf.ndim != 2 || buf.shape[1] != 2) && (buf.ndim != 1 || buf.shape[0] != 0))
    throw std::runtime_error("Diagram must be an array of size n x 2");
  // In the case of shape (0), avoid reading non-existing strides[1] even if we won't use it.
  py::ssize_t stride1 = buf.ndim == 2 ? buf.strides[1] : 0;
  auto cnt = boost::counting_range<py::ssize_t>(0, buf.shape[0]);

  char* p = static_cast<char*>(buf.ptr);
  auto h = buf.strides[0];
  auto w = stride1;
  // Get m[i,0] and m[i,1] as a pair
  auto pairify = [=](py::ssize_t i) {
    char* birth = p + i * h;
    char* death = birth + w;
    return build_point(*(double*)birth, *(double*)death, i);
  };
  return boost::adaptors::transform(cnt, pairify);
  // Be careful that the returned range cannot contain references to dead temporaries.
}
