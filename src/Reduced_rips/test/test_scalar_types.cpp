/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE reduced_rips_scalar_types
#include <boost/test/unit_test.hpp>

#include <array>
#include <type_traits>
#include <vector>

#include "reduced_rips_test_fixtures.h"

// The reduction is templated on Filtration_value: distances, edge lengths, lune deaths and the barcode are all
// carried in that type. A given instantiation must recover the same degree-1 barcode as the double ground
// truth, to that type's precision, on both the point-cloud and the distance-matrix paths.
template <class FV>
void check_scalar_type(const char* tag, double tol) {
  using RR = Gudhi::reduced_rips::Reduced_rips<FV>;
  static_assert(std::is_same_v<typename RR::Filtration_value, FV>, "Reduced_rips<FV>::Filtration_value must be FV");
  static_assert(std::is_same_v<typename RR::Persistence_interval, std::array<FV, 2>>,
                "the barcode must be arrays of FV, with no widening to double");
  BOOST_TEST_CONTEXT(tag) {
    auto pts = circle(40, 2);
    Bars truth = full_rips_h1(pts);

    std::vector<std::array<FV, 2>> from_points = RR::from_points(pts).persistence();  // copy out of the temporary
    std::vector<std::array<FV, 2>> from_matrix = RR::from_distance_matrix(full_distance_matrix(pts)).persistence();

    BOOST_CHECK(bars_close(truth, to_double_bars(from_points), tol));
    BOOST_CHECK(bars_close(truth, to_double_bars(from_matrix), tol));
  }
}

BOOST_AUTO_TEST_CASE(filtration_value_scalar_types) {
  // float carries ~7 significant digits (squared then square-rooted), so it needs a looser tolerance than the
  // wider types; all three must still recover the loop.
  check_scalar_type<float>("float", 1e-3);
  check_scalar_type<double>("double", 1e-6);
  check_scalar_type<long double>("long double", 1e-9);
}

BOOST_AUTO_TEST_CASE(integer_distance_matrix) {
  // The matrix geometry carries the supplied dissimilarities verbatim (no squaring, no square root), so an
  // integer Filtration_value keeps the whole reduction in exact integer arithmetic. A 4-cycle with edge length
  // 2 and diagonals 3 (a scaled square) has a single H1 loop: born when the four length-2 edges close the cycle,
  // dying when a diagonal triangle fills it, giving the exact bar (2, 3).
  using RR = Gudhi::reduced_rips::Reduced_rips<int>;
  static_assert(std::is_same_v<RR::Filtration_value, int>);
  std::vector<std::vector<int>> lower = {
      {},         // point 0
      {2},        // d(1,0)
      {3, 2},     // d(2,0), d(2,1)
      {2, 3, 2},  // d(3,0), d(3,1), d(3,2)
  };
  std::vector<std::array<int, 2>> bc = RR::from_distance_matrix(lower).persistence();
  BOOST_REQUIRE_EQUAL(bc.size(), 1U);
  BOOST_CHECK_EQUAL(bc.front()[0], 2);  // exact integer birth
  BOOST_CHECK_EQUAL(bc.front()[1], 3);  // exact integer death
}
