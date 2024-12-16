/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

namespace py = nanobind;
typedef py::ndarray<double> Dgm;

// build_point(double birth, double death, ssize_t index) -> Point
template<class BuildPoint>
inline auto numpy_to_range_of_pairs(Dgm dgm, BuildPoint build_point) {
  // shape (n,2) or (0) for empty
  if((dgm.ndim()!=2 || dgm.shape(1)!=2) && (dgm.ndim()!=1 || dgm.shape(0)!=0))
    throw std::runtime_error("Diagram must be an array of size n x 2");
  // In the case of shape (0), avoid reading non-existing strides[1] even if we won't use it.
  py::ssize_t stride1 = dgm.ndim() == 2 ? dgm.stride(1) : 0;
  auto cnt = boost::counting_range<py::ssize_t>(0, dgm.shape(0));

  auto ptr = dgm.shape_ptr();
  auto h = dgm.stride(0);
  auto w = stride1;
  // Get m[i,0] and m[i,1] as a pair
  auto pairify = [=](py::ssize_t i){
    auto birth = ptr + i * h;
    auto death = birth + w;
    return build_point(*(double*)birth, *(double*)death, i);
  };
  return boost::adaptors::transform(cnt, pairify);
  // Be careful that the returned range cannot contain references to dead temporaries.
}
