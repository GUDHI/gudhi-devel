/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2025/03 Vincent Rouvreau & Hannah Schreiber: Use nanobind instead of Cython for python bindings and renaming.
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstddef>  //std::size_t

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

using Tensor_dgm = nanobind::ndarray<const double, nanobind::shape<-1, 2> >;
using Tensor_dgm_view = decltype(std::declval<Tensor_dgm>().view());

// build_point(double birth, double death, size_t index) -> Point
template <class BuildPoint>
inline auto array_to_range_of_pairs(const Tensor_dgm_view& dgm, BuildPoint&& build_point)
{
  auto cnt = boost::counting_range<std::size_t>(0, dgm.shape(0));
  // Get m[i,0] and m[i,1] as a pair
  auto pairify = [&](nanobind::ssize_t i) {
    double birth = dgm(i, 0);
    double death = dgm(i, 1);
    return build_point(birth, death, i);
  };
  return boost::adaptors::transform(cnt, pairify);
  // Be careful that the returned range cannot contain references to dead temporaries.
}

