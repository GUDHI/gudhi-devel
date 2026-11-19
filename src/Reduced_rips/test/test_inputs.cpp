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
#define BOOST_TEST_MODULE reduced_rips_inputs
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

#include "reduced_rips_test_fixtures.h"

// ---- Distance-matrix backend ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(distance_matrix_matches_point_cloud) {
  // The distance-matrix path, fed a Euclidean cloud's own distances, must reproduce the coordinate path.
  for (unsigned dim : {2U, 3U, 5U}) {
    BOOST_TEST_CONTEXT("dim=" << dim) {
      auto pts = circle(50, dim);
      auto from_points = Reduced_rips::from_points(pts);
      Reduced_rips from_matrix = Reduced_rips::from_distance_matrix(full_distance_matrix(pts));
      BOOST_REQUIRE(!from_points.persistence().empty());  // the circle has H1, so equal-but-empty is a failure
      BOOST_CHECK(bars_close(from_points.persistence(), from_matrix.persistence()));
    }
  }
}

BOOST_AUTO_TEST_CASE(full_and_lower_triangular_agree) {
  // The same square given as a full symmetric matrix yields the same barcode as the lower-triangular form.
  const double s2 = std::sqrt(2.0);
  std::vector<std::vector<double>> lower = {{}, {1.0}, {s2, 1.0}, {1.0, s2, 1.0}};
  std::vector<std::vector<double>> full = {
      {0.0, 1.0, s2, 1.0}, {1.0, 0.0, 1.0, s2}, {s2, 1.0, 0.0, 1.0}, {1.0, s2, 1.0, 0.0}};
  auto a = Reduced_rips::from_distance_matrix(lower).persistence();
  auto b = Reduced_rips::from_distance_matrix(full).persistence();
  BOOST_REQUIRE(!a.empty());  // the square has one bar, so equal-but-empty is a failure
  BOOST_CHECK(bars_close(a, b));
}

BOOST_AUTO_TEST_CASE(non_metric_dissimilarity_is_accepted) {
  // The reduction never uses the triangle inequality, so a symmetric matrix that badly violates it is still
  // a valid input and must compute the exact degree-1 barcode of that dissimilarity's Vietoris-Rips
  // filtration. Compare against the full-complex ground truth built from the same matrix.
  std::vector<std::vector<double>> full = {
      {0.0, 5.0, 1.0, 1.0}, {5.0, 0.0, 1.0, 1.0}, {1.0, 1.0, 0.0, 5.0}, {1.0, 1.0, 5.0, 0.0}};
  Bars truth = full_rips_h1_from_matrix(full);
  Bars reduced = Reduced_rips::from_distance_matrix(full).persistence();
  BOOST_CHECK(bars_close(truth, reduced));
}

// ---- Search-strategy and neighbor-budget invariance ----------------------------------------------------

BOOST_AUTO_TEST_CASE(search_strategies_agree) {
  // kd-tree and brute-force neighbor search are an implementation choice; the barcode must not depend on it,
  // in low ambient dimension (where automatic picks kd-tree) and high (where it picks brute-force).
  for (std::size_t dim : {std::size_t{3}, std::size_t{6}}) {
    BOOST_TEST_CONTEXT("dim=" << dim) {
      auto pts = random_cloud(50, dim, 11);
      Bars various = Reduced_rips::from_points(pts, 0, Reduced_rips::Search::automatic).persistence();
      Bars kd = Reduced_rips::from_points(pts, 0, Reduced_rips::Search::kd_tree).persistence();
      Bars brute = Reduced_rips::from_points(pts, 0, Reduced_rips::Search::brute_force).persistence();
      BOOST_REQUIRE(!various.empty());  // these fixed random clouds have H1; equal-but-empty is a failure
      BOOST_CHECK(bars_close(various, kd));
      BOOST_CHECK(bars_close(various, brute));
    }
  }
}

BOOST_AUTO_TEST_CASE(initial_neighbor_budget_does_not_change_result) {
  // num_neighbors only seeds the heap; the frontier grows on demand, so the barcode is independent of it.
  auto pts = random_cloud(50, 3, 21);
  Bars budget_default = Reduced_rips::from_points(pts, 0).persistence();
  Bars budget_small = Reduced_rips::from_points(pts, 3).persistence();
  Bars budget_large = Reduced_rips::from_points(pts, 40).persistence();
  BOOST_REQUIRE(!budget_default.empty());  // this fixed random cloud has H1; equal-but-empty is a failure
  BOOST_CHECK(bars_close(budget_default, budget_small));
  BOOST_CHECK(bars_close(budget_default, budget_large));
}

#ifdef GUDHI_USE_TBB
#include <tbb/global_control.h>

BOOST_AUTO_TEST_CASE(single_thread_matches_parallel) {
  // The lune evaluations are parallel but the reduction consumes their results in pop order, so the barcode
  // must be independent of the number of TBB workers. Pin TBB to one worker and compare with the default run.
  auto pts = torus(12, 8);
  Bars parallel = Reduced_rips::from_points(pts).persistence();
  Bars serial;
  {
    tbb::global_control gc(tbb::global_control::max_allowed_parallelism, 1);
    serial = Reduced_rips::from_points(pts).persistence();
  }
  BOOST_REQUIRE(!parallel.empty());
  BOOST_CHECK(bars_close(parallel, serial));
}
#endif
