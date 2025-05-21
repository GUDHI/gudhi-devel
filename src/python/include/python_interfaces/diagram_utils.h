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

#include <cstddef>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

using Tensor_dgm = nanobind::ndarray<const double, nanobind::shape<-1, 2> >;
using Sequence_dgm = nanobind::sequence;
using List_dgm = nanobind::list;

// build_point(double birth, double death, size_t index) -> Point
template <class BuildPoint>
inline auto array_to_range_of_pairs(const Tensor_dgm& dgm, BuildPoint&& build_point)
{
  auto cnt = boost::counting_range<std::size_t>(0, dgm.shape(0));
  // Get m[i,0] and m[i,1] as a pair
  auto pairify = [=](nanobind::ssize_t i) {
    double birth = dgm(i, 0);
    double death = dgm(i, 1);
    return build_point(birth, death, i);
  };
  return boost::adaptors::transform(cnt, pairify);
  // Be careful that the returned range cannot contain references to dead temporaries.
}

template <class Array, class BuildPoint>
inline auto _array_to_range_of_pairs(const Array& dgm, BuildPoint&& build_point)
{
  auto size = nanobind::len(dgm);
  std::vector<decltype(build_point(0, 0, 0))> cnt(size);
  std::size_t i = 0;
  for (auto it = dgm.begin(); it != dgm.end(); ++it) {
    const auto& p = *it;
    auto itP = p.begin();
    if (itP == p.end()) throw std::runtime_error("Diagram must be an array of size n x 2");
    double birth = nanobind::cast<double>(*itP);
    ++itP;
    if (itP == p.end()) throw std::runtime_error("Diagram must be an array of size n x 2");
    double death = nanobind::cast<double>(*itP);
    ++itP;
    if (itP != p.end()) throw std::runtime_error("Diagram must be an array of size n x 2");
    cnt[i] = build_point(birth, death, i);
    ++i;
  }
  return cnt;
}

// build_point(double birth, double death, size_t index) -> Point
template <class BuildPoint>
inline auto array_to_range_of_pairs(const Sequence_dgm& dgm, BuildPoint&& build_point)
{
  return _array_to_range_of_pairs(dgm, std::forward<BuildPoint>(build_point));
}

// build_point(double birth, double death, size_t index) -> Point
template <class BuildPoint>
inline auto array_to_range_of_pairs(const List_dgm& dgm, BuildPoint&& build_point)
{
  return _array_to_range_of_pairs(dgm, std::forward<BuildPoint>(build_point));
}
