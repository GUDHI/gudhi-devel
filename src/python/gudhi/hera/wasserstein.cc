/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <pybind11_diagram_utils.h>

#ifdef _MSC_VER
// https://github.com/grey-narn/hera/issues/3
// ssize_t is a non-standard type (well, posix)
using py::ssize_t;
#endif

#include <hera/wasserstein.h>
#include <gudhi/Debug_utils.h>

// Unlike bottleneck, for wasserstein, we need to add the index ourselves (if we want the matching)
static auto make_hera_point(double x, double y, py::ssize_t i) { return hera::DiagramPoint<double>(x, y, i); };

py::object wasserstein_distance(
    Dgm d1, Dgm d2,
    double wasserstein_power, double internal_p,
    double delta, bool return_matching)
{
  // I *think* the call to request() in numpy_to_range_of_pairs has to be before releasing the GIL.
  auto diag1 = numpy_to_range_of_pairs(d1, make_hera_point);
  auto diag2 = numpy_to_range_of_pairs(d2, make_hera_point);
  int n1 = boost::size(diag1);
  int n2 = boost::size(diag2);
  hera::AuctionResult<double> res;
  double dist;

  { // No Python allowed in this section
    py::gil_scoped_release release;

    hera::AuctionParams<double> params;
    params.wasserstein_power = wasserstein_power;
    // hera encodes infinity as -1...
    if(std::isinf(internal_p)) internal_p = hera::get_infinity<double>();
    params.internal_p = internal_p;
    params.delta = delta;
    if(return_matching) {
      params.return_matching = true;
      params.match_inf_points = true;
    }
    // The extra parameters are purposely not exposed for now.
    res = hera::wasserstein_cost_detailed(diag1, diag2, params);
    dist = std::pow(res.cost, 1./params.wasserstein_power);
  }

  if(!return_matching)
    return py::cast(dist);

  if(dist == std::numeric_limits<double>::infinity())
    return py::make_tuple(dist, py::none());

  // bug in Hera, matching_a_to_b_ is empty if one diagram is empty or both diagrams contain the same points
  if(res.matching_a_to_b_.size() == 0) {
    if(n1 == 0) { // diag1 is empty
      py::array_t<int> matching({{ n2, 2 }}, nullptr);
      auto m = matching.mutable_unchecked();
      for(int j=0; j<n2; ++j){
        m(j, 0) = -1;
        m(j, 1) = j;
      }
      return py::make_tuple(dist, matching);
    }
    if(n2 == 0) { // diag2 is empty
      py::array_t<int> matching({{ n1, 2 }}, nullptr);
      auto m = matching.mutable_unchecked();
      for(int i=0; i<n1; ++i){
        m(i, 0) = i;
        m(i, 1) = -1;
      }
      return py::make_tuple(dist, matching);
    }
    // The only remaining case should be that the 2 diagrams are identical, but possibly shuffled
    GUDHI_CHECK(n1==n2, "unexpected bug in Hera?");
    std::vector v1(boost::begin(diag1), boost::end(diag1));
    std::vector v2(boost::begin(diag2), boost::end(diag2));
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    py::array_t<int> matching({{ n1, 2 }}, nullptr);
    auto m = matching.mutable_unchecked();
    for(int i=0; i<n1; ++i){
      GUDHI_CHECK(v1[i][0]==v2[i][0] && v1[i][1]==v2[i][1], "unexpected bug in Hera?");
      m(i, 0) = v1[i].get_id();
      m(i, 1) = v2[i].get_id();
    }
    return py::make_tuple(dist, matching);

  }

  // bug in Hera, diagonal points are ignored and don't appear in matching_a_to_b_
  for(auto p : diag1)
    if(p[0] == p[1]) { auto id = p.get_id(); res.matching_a_to_b_[id] = -id-1; }
  for(auto p : diag2)
    if(p[0] == p[1]) { auto id = p.get_id(); res.matching_a_to_b_[-id-1] = id; }

  py::array_t<int> matching({{ n1 + n2, 2 }}, nullptr);
  auto m = matching.mutable_unchecked();
  int cur = 0;
  for(auto x : res.matching_a_to_b_){
    if(x.first < 0) {
      if(x.second < 0) {
      } else {
        m(cur, 0) = -1;
        m(cur, 1) = x.second;
        ++cur;
      }
    } else {
      if(x.second < 0) {
        m(cur, 0) = x.first;
        m(cur, 1) = -1;
        ++cur;
      } else {
        m(cur, 0) = x.first;
        m(cur, 1) = x.second;
        ++cur;
      }
    }
  }
  // n1+n2 was too much, it only happens if everything matches to the diagonal, so we return matching[:cur,:]
  py::array_t<int> ret({{ cur, 2 }}, {{ matching.strides()[0], matching.strides()[1] }}, matching.data(), matching);
  return py::make_tuple(dist, ret);
}

PYBIND11_MODULE(wasserstein, m) {
      m.def("wasserstein_distance", &wasserstein_distance,
          py::arg("X"), py::arg("Y"),
          py::arg("order") = 1,
          py::arg("internal_p") = std::numeric_limits<double>::infinity(),
          py::arg("delta") = .01,
          py::arg("matching") = false,
          R"pbdoc(
        Compute the Wasserstein distance between two diagrams.
        Points at infinity are supported.

        Parameters:
            X (n x 2 numpy array): First diagram
            Y (n x 2 numpy array): Second diagram
            order (float): Wasserstein exponent W_q
            internal_p (float): Internal Minkowski norm L^p in R^2
            delta (float): Relative error 1+delta
            matching (bool): if ``True``, computes and returns the optimal matching between X and Y, encoded as a (n x 2) np.array [...[i,j]...], meaning the i-th point in X is matched to the j-th point in Y, with the convention that (-1) represents the diagonal. If the distance between two diagrams is +inf (which happens if the cardinalities of essential parts differ) and the matching is requested, it will be set to ``None`` (any matching is optimal).

        Returns:
            float|Tuple[float,numpy.array|None]: Approximate Wasserstein distance W_q(X,Y), and optionally the corresponding matching
    )pbdoc");
}
