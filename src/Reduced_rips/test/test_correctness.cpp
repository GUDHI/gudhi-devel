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
#define BOOST_TEST_MODULE reduced_rips_correctness
#include <boost/test/unit_test.hpp>

#include <gudhi/Unitary_tests_utils.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "reduced_rips_test_fixtures.h"

// ---- Correctness against the full Vietoris-Rips filtration ---------------------------------------------

BOOST_AUTO_TEST_CASE(matches_full_vietoris_rips_h1) {
  // The headline guarantee of the reference paper: the reduced filtration has the same degree-1 barcode as
  // the full Vietoris-Rips filtration. Check it exactly on a spread of inputs, exercising the 2D and 3D
  // Delaunay relative-neighborhood-graph paths and the dimension-free O(n^2) path, and both backends.
  struct Case {
    const char* name;
    Cloud pts;
  };
  std::vector<Case> cases = {
      {"circle 2D", circle(40, 2)},
      {"circle 5D", circle(36, 5)},
      {"torus 3D", torus(12, 8)},
      {"two loops 2D", two_circles(30)},
      {"random 2D", random_cloud(60, 2, 1)},
      {"random 3D", random_cloud(55, 3, 2)},
      {"random 5D", random_cloud(45, 5, 3)},
  };
  for (const auto& c : cases) {
    BOOST_TEST_CONTEXT(c.name) {
      Bars truth = full_rips_h1(c.pts);
      Bars from_points = Reduced_rips::from_points(c.pts).persistence();
      Bars from_matrix = Reduced_rips::from_distance_matrix(full_distance_matrix(c.pts)).persistence();
      BOOST_CHECK(bars_close(truth, from_points));
      BOOST_CHECK(bars_close(truth, from_matrix));
    }
  }
}

// ---- Topology recovery on hand-understood inputs -------------------------------------------------------

BOOST_AUTO_TEST_CASE(circle_has_one_dominant_loop) {
  for (std::size_t dim : {std::size_t{2}, std::size_t{5}}) {  // 2D Delaunay path and >=4D general path
    BOOST_TEST_CONTEXT("dim=" << dim) {
      Bars bc = Reduced_rips::from_points(circle(60, dim)).persistence();
      // Exactly one prominent loop, and it dies past 1.5 (near the unit circle's diameter of 2).
      BOOST_REQUIRE(!bc.empty());  // max_element below dereferences the range
      BOOST_CHECK_EQUAL(count_prominent(bc, 0.5), 1u);
      auto loop = *std::max_element(bc.begin(), bc.end(), [](const std::array<double, 2>& a,
                                                             const std::array<double, 2>& b) {
        return a[1] - a[0] < b[1] - b[0];
      });
      BOOST_CHECK_GT(loop[1], 1.5);
    }
  }
}

BOOST_AUTO_TEST_CASE(two_separated_circles_give_two_loops) {
  Bars bc = Reduced_rips::from_points(two_circles(30)).persistence();
  BOOST_CHECK_EQUAL(count_prominent(bc, 0.5), 2u);
}

BOOST_AUTO_TEST_CASE(unit_square_known_bar) {
  // Four points of a unit square (sides 1, diagonals sqrt 2): a single H1 loop born when the four unit edges
  // close the cycle and dying when a diagonal fills it in.
  // Use the matrix form to pin exact distances independent of trig rounding.
  const double s2 = std::sqrt(2.0);
  std::vector<std::vector<double>> lower = {
      {},             // point 0
      {1.0},          // d(1,0)
      {s2, 1.0},      // d(2,0), d(2,1)
      {1.0, s2, 1.0}  // d(3,0), d(3,1), d(3,2)
  };
  Bars bc = Reduced_rips::from_distance_matrix(lower).persistence();
  BOOST_REQUIRE_EQUAL(bc.size(), 1U);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(bc.front()[0], 1.0, 1e-6);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(bc.front()[1], s2, 1e-6);
}
