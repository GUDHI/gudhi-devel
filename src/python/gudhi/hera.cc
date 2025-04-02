/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2025/01 Vincent Rouvreau: Use nanobind instead of PyBind11 for python bindings. Merge of bottleneck.cc and
 *          wasserstein.cc, because nanobind doesn't accept 2 modules named hera
 *      - YYYY/MM Author: Description of the modification
 */

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <hera/bottleneck.h>
#include <hera/wasserstein.h>

#include <gudhi/Debug_utils.h>
#include <python_interfaces/diagram_utils.h>

namespace nb = nanobind;

#ifdef _MSC_VER
// https://github.com/grey-narn/hera/issues/3
// ssize_t is a non-standard type (well, posix)
using nb::ssize_t;
#endif

// Indices are added internally in bottleneck_distance, they are not needed in the input.
static auto make_point(double x, double y, nb::ssize_t) { return std::pair(x, y); };

double bottleneck_distance(const Dgm& d1, const Dgm& d2, double delta)
{
  // I *think* the call to request() in array_to_range_of_pairs has to be before releasing the GIL.
  auto diag1 = array_to_range_of_pairs(d1, make_point);
  auto diag2 = array_to_range_of_pairs(d2, make_point);

  nb::gil_scoped_release release;

  if (delta == 0)
    return hera::bottleneckDistExact(diag1, diag2);
  else
    return hera::bottleneckDistApprox(diag1, diag2, delta);
}

// Unlike bottleneck, for wasserstein, we need to add the index ourselves (if we want the matching)
static auto make_hera_point(double x, double y, nb::ssize_t i) { return hera::DiagramPoint<double>(x, y, i); };

nb::object wasserstein_distance(const Dgm& d1,
                                const Dgm& d2,
                                double wasserstein_power,
                                double internal_p,
                                double delta,
                                bool return_matching)
{
  auto diag1 = array_to_range_of_pairs(d1, make_hera_point);
  auto diag2 = array_to_range_of_pairs(d2, make_hera_point);
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
      std::vector<int> arr(2 * n2);
      nb::ndarray<int, nb::numpy, nb::ndim<2>> matching(arr.data(), {static_cast<size_t>(n2), 2});
      for (int j = 0; j < n2; ++j) {
        matching(j, 0) = -1;
        matching(j, 1) = j;
      }
      return nb::make_tuple(dist, matching);
    }
    if (n2 == 0) {  // diag2 is empty
      std::vector<int> arr(2 * n1);
      nb::ndarray<int, nb::numpy, nb::ndim<2>> matching(arr.data(), {static_cast<size_t>(n1), 2});
      for (int i = 0; i < n1; ++i) {
        matching(i, 0) = i;
        matching(i, 1) = -1;
      }
      return nb::make_tuple(dist, matching);
    }
    // The only remaining case should be that the 2 diagrams are identical, but possibly shuffled
    GUDHI_CHECK(n1 == n2, "unexpected bug in Hera?");
    std::vector v1(boost::begin(diag1), boost::end(diag1));
    std::vector v2(boost::begin(diag2), boost::end(diag2));
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::vector<int> arr(2 * n1);
    nb::ndarray<int, nb::numpy, nb::ndim<2>> matching(arr.data(), {static_cast<size_t>(n1), 2});
    for (int i = 0; i < n1; ++i) {
      GUDHI_CHECK(v1[i][0] == v2[i][0] && v1[i][1] == v2[i][1], "unexpected bug in Hera?");
      matching(i, 0) = v1[i].get_id();
      matching(i, 1) = v2[i].get_id();
    }
    return nb::make_tuple(dist, matching);
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

  std::vector<int> arr(2 * (n1 + n2));
  int cur = 0;
  for (auto x : res.matching_a_to_b_) {
    if (x.first < 0) {
      if (x.second < 0) {
      } else {
        arr[2 * cur + 0] = -1;
        arr[2 * cur + 1] = x.second;
        cur++;
      }
    } else {
      if (x.second < 0) {
        arr[2 * cur + 0] = x.first;
        arr[2 * cur + 1] = -1;
        cur++;
      } else {
        arr[2 * cur + 0] = x.first;
        arr[2 * cur + 1] = x.second;
        cur++;
      }
    }
  }
  // n1+n2 was too much, it only happens if everything matches to the diagonal, so we return matching[:cur,:]
  arr.resize(2 * cur);
  nb::ndarray<int, nb::numpy, nb::ndim<2>> ret(arr.data(), {static_cast<size_t>(cur), 2});
  return nb::make_tuple(dist, ret);
}

NB_MODULE(hera, m)
{
  m.def("bottleneck_distance",
        &bottleneck_distance,
        nb::arg("X"),
        nb::arg("Y"),
        nb::arg("delta") = .01,
        R"pbdoc(
    Compute the Bottleneck distance between two diagrams.
    Points at infinity are supported.

    .. note::
       Points on the diagonal are not supported and must be filtered out before calling this function.

    Parameters:
        X (n x 2 numpy array): First diagram
        Y (n x 2 numpy array): Second diagram
        delta (float): Relative error 1+delta

    Returns:
        float: (approximate) bottleneck distance d_B(X,Y)
  )pbdoc");
  m.def("wasserstein_distance",
        &wasserstein_distance,
        nb::arg("X"),
        nb::arg("Y"),
        nb::arg("order") = 1,
        nb::arg("internal_p") = std::numeric_limits<double>::infinity(),
        nb::arg("delta") = .01,
        nb::arg("matching") = false,
        R"pbdoc(
    Compute the Wasserstein distance between two diagrams.
    Points at infinity are supported.

    Parameters:
        X (n x 2 numpy array): First diagram
        Y (n x 2 numpy array): Second diagram
        order (float): Wasserstein exponent W_q
        internal_p (float): Internal Minkowski norm L^p in R^2
        delta (float): Relative error 1+delta
        matching (bool): if ``True``, computes and returns the optimal matching between X and Y, encoded as a
            (n x 2) np.array [...[i,j]...], meaning the i-th point in X is matched to the j-th point in Y, with the
            convention that (-1) represents the diagonal. If the distance between two diagrams is +inf (which happens
            if the cardinalities of essential parts differ) and the matching is requested, it will be set to ``None``
            (any matching is optimal).

        Returns:
            float|Tuple[float,numpy.array|None]: Approximate Wasserstein distance W_q(X,Y), and optionally the
                corresponding matching
  )pbdoc");
}
