/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2025/01 Vincent Rouvreau: Use nanobind instead of PyBind11 for python bindings. Merge of bottleneck.cc and
 *          wasserstein.cc, because nanobind doesn't accept 2 modules named hera
 *      - 2025/04 Hannah Schreiber: Re-add possibility of native Python sequences as input
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstdint>  // missing in hera/wasserstein.h

#include <nanobind/nanobind.h>

#ifdef _MSC_VER
// https://github.com/grey-narn/hera/issues/3
// ssize_t is a non-standard type (well, posix)
using nanobind::ssize_t;
#endif

#include <hera/bottleneck.h>
#include <hera/wasserstein.h>

#include <gudhi/Debug_utils.h>
#include <python_interfaces/diagram_utils.h>
#include <python_interfaces/numpy_utils.h>

namespace nb = nanobind;

// Indices are added internally in bottleneck_distance, they are not needed in the input.
static auto _make_point(double x, double y, std::size_t) { return std::pair(x, y); };

template <class Dgm>
double _bottleneck_distance(const Dgm& d1, const Dgm& d2, double delta)
{
  // I *think* the call to request() in array_to_range_of_pairs has to be before releasing the GIL.
  auto diag1 = array_to_range_of_pairs(d1, _make_point);
  auto diag2 = array_to_range_of_pairs(d2, _make_point);

  nb::gil_scoped_release release;

  if (delta == 0)
    return hera::bottleneckDistExact(diag1, diag2);
  else
    return hera::bottleneckDistApprox(diag1, diag2, delta);
}

// Unlike bottleneck, for wasserstein, we need to add the index ourselves (if we want the matching)
static auto _make_hera_point(double x, double y, std::size_t i) { return hera::DiagramPoint<double>(x, y, i); };

template <class Dgm>
nb::object _wasserstein_distance(const Dgm& d1,
                                 const Dgm& d2,
                                 double wasserstein_power,
                                 double internal_p,
                                 double delta,
                                 bool return_matching)
{
  auto diag1 = array_to_range_of_pairs(d1, _make_hera_point);
  auto diag2 = array_to_range_of_pairs(d2, _make_hera_point);
  int n1 = boost::size(diag1);
  int n2 = boost::size(diag2);
  hera::AuctionResult<double> res;
  double dist;

  {  // No Python allowed in this section
    nb::gil_scoped_release release;

    hera::AuctionParams<double> params;
    params.wasserstein_power = wasserstein_power;
    // hera encodes infinity as -1...
    if (std::isinf(internal_p)) internal_p = hera::get_infinity<double>();
    params.internal_p = internal_p;
    params.delta = delta;
    if (return_matching) {
      params.return_matching = true;
      params.match_inf_points = true;
    }
    // The extra parameters are purposely not exposed for now.
    res = hera::wasserstein_cost_detailed(diag1, diag2, params);
    dist = std::pow(res.cost, 1. / params.wasserstein_power);
  }

  if (!return_matching) return nb::cast(dist);

  if (dist == std::numeric_limits<double>::infinity()) return nb::make_tuple(dist, nb::none());

  // bug in Hera, matching_a_to_b_ is empty if one diagram is empty or both diagrams contain the same points
  if (res.matching_a_to_b_.size() == 0) {
    if (n1 == 0) {  // diag1 is empty
      auto matching = new int[2 * n2];
      int i = 0;
      for (int j = 0; j < n2; ++j) {
        matching[i] = -1;
        matching[i + 1] = j;
        i += 2;
      }
      return nb::make_tuple(dist, _wrap_as_numpy_array(matching, n2, 2));
    }
    if (n2 == 0) {  // diag2 is empty
      auto matching = new int[2 * n1];
      int i = 0;
      for (int j = 0; j < n1; ++j) {
        matching[i] = j;
        matching[i + 1] = -1;
        i += 2;
      }
      return nb::make_tuple(dist, _wrap_as_numpy_array(matching, n1, 2));
    }
    // The only remaining case should be that the 2 diagrams are identical, but possibly shuffled
    GUDHI_CHECK(n1 == n2, "unexpected bug in Hera?");
    std::vector v1(boost::begin(diag1), boost::end(diag1));
    std::vector v2(boost::begin(diag2), boost::end(diag2));
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    auto matching = new int[2 * n1];
    int j = 0;
    for (int i = 0; i < n1; ++i) {
      GUDHI_CHECK(v1[i][0] == v2[i][0] && v1[i][1] == v2[i][1], "unexpected bug in Hera?");
      matching[j] = v1[i].get_id();
      matching[j + 1] = v2[i].get_id();
      j += 2;
    }
    return nb::make_tuple(dist, _wrap_as_numpy_array(matching, n1, 2));
  }

  // bug in Hera, diagonal points are ignored and don't appear in matching_a_to_b_
  for (auto p : diag1)
    if (p[0] == p[1]) {
      auto id = p.get_id();
      res.matching_a_to_b_[id] = -id - 1;
    }
  for (auto p : diag2)
    if (p[0] == p[1]) {
      auto id = p.get_id();
      res.matching_a_to_b_[-id - 1] = id;
    }

  auto matching = new int[2 * (n1 + n2)];
  int cur = 0;
  for (auto x : res.matching_a_to_b_) {
    if (x.first < 0) {
      if (x.second < 0) {
      } else {
        matching[2 * cur + 0] = -1;
        matching[2 * cur + 1] = x.second;
        cur++;
      }
    } else {
      if (x.second < 0) {
        matching[2 * cur + 0] = x.first;
        matching[2 * cur + 1] = -1;
        cur++;
      } else {
        matching[2 * cur + 0] = x.first;
        matching[2 * cur + 1] = x.second;
        cur++;
      }
    }
  }
  // n1+n2 was too much, it only happens if everything matches to the diagonal, so we return matching[:cur,:]
  return nb::make_tuple(dist, _wrap_as_numpy_array(matching, cur, 2));
}

NB_MODULE(_hera_ext, m)
{
  m.attr("__license__") = "BSD 3-Clause";
  m.def("_bottleneck_distance_tensor",
        &_bottleneck_distance<Tensor_dgm>,
        nb::arg("X"),
        nb::arg("Y"),
        nb::arg("delta") = .01);
  m.def(
      "_bottleneck_distance_list", &_bottleneck_distance<List_dgm>, nb::arg("X"), nb::arg("Y"), nb::arg("delta") = .01);
  m.def("_bottleneck_distance_sequence",
        &_bottleneck_distance<Sequence_dgm>,
        nb::arg("X"),
        nb::arg("Y"),
        nb::arg("delta") = .01);
  m.def("_wasserstein_distance_tensor",
        &_wasserstein_distance<Tensor_dgm>,
        nb::arg("X"),
        nb::arg("Y"),
        nb::arg("order") = 1,
        nb::arg("internal_p") = std::numeric_limits<double>::infinity(),
        nb::arg("delta") = .01,
        nb::arg("matching") = false);
  m.def("_wasserstein_distance_list",
        &_wasserstein_distance<List_dgm>,
        nb::arg("X"),
        nb::arg("Y"),
        nb::arg("order") = 1,
        nb::arg("internal_p") = std::numeric_limits<double>::infinity(),
        nb::arg("delta") = .01,
        nb::arg("matching") = false);
  m.def("_wasserstein_distance_sequence",
        &_wasserstein_distance<Sequence_dgm>,
        nb::arg("X"),
        nb::arg("Y"),
        nb::arg("order") = 1,
        nb::arg("internal_p") = std::numeric_limits<double>::infinity(),
        nb::arg("delta") = .01,
        nb::arg("matching") = false);
}
