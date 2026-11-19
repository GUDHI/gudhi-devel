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
#define BOOST_TEST_MODULE reduced_rips_edge_cases
#include <boost/test/unit_test.hpp>

#include <gudhi/Points_off_io.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "reduced_rips_test_fixtures.h"

// ---- Points exactly on lune boundaries -----------------------------------------------------------------

BOOST_AUTO_TEST_CASE(equilateral_triangle_boundary_is_trivial) {
  // Three points at mutual distance 1: each vertex lies *exactly* on the lune boundary of the opposite edge
  // (the d(a,c) == d(b,c) == d(a,b) double-boundary case). The true relative neighborhood graph uses the open
  // lune, so a boundary point does not remove an edge: all three edges are kept, the single cycle is filled by
  // the triangle at the same filtration value, and H1 is trivial. The matrix pins the distances bit-exact, so
  // the boundary comparisons are genuine equalities rather than near-ties.
  std::vector<std::vector<double>> full = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
  Bars reduced = Reduced_rips::from_distance_matrix(full).persistence();
  BOOST_CHECK_EQUAL(finite_positive(reduced).size(), 0u);
  BOOST_CHECK(bars_close(full_rips_h1_from_matrix(full), reduced));
}

BOOST_AUTO_TEST_CASE(boundary_point_on_loop_matches_full_rips) {
  // A unit-square loop (bar (1, sqrt2)) plus a fifth point sitting *exactly* on a lune boundary of one square
  // edge: distance 1 (== the edge length) from corner 0 and strictly less (0.5) from corner 1. This exercises
  // in_lune's single-boundary branch in the reduction while the true RNG (open lune) still keeps the edge, so
  // the loop must survive. Distances are pinned exactly by the matrix; the reduced barcode must match the full
  // Vietoris-Rips ground truth on the same matrix.
  const double s2 = std::sqrt(2.0);
  std::vector<std::vector<double>> full = {{0.0, 1.0, s2, 1.0, 1.0},
                                           {1.0, 0.0, 1.0, s2, 0.5},
                                           {s2, 1.0, 0.0, 1.0, 2.0},
                                           {1.0, s2, 1.0, 0.0, 2.0},
                                           {1.0, 0.5, 2.0, 2.0, 0.0}};
  Bars truth = full_rips_h1_from_matrix(full);
  Bars reduced = Reduced_rips::from_distance_matrix(full).persistence();
  BOOST_CHECK(bars_close(truth, reduced));
  BOOST_CHECK_EQUAL(count_prominent(reduced, 0.3), count_prominent(truth, 0.3));
}

BOOST_AUTO_TEST_CASE(rng_cycle_rank_paths_agree_on_boundary_ties) {
  // A 3-4-5 configuration whose apex lies bit-exactly on a lune boundary: d(0,1) == d(0,2) == 5 with
  // d(1,2) < 5, all representable in double. The strict open lune keeps a boundary-tied edge, so the true
  // RNG has all three edges and cycle rank 1, on every path. The general/matrix supergraph build once broke
  // the (0,1)/(0,2) tie by index and silently dropped an edge (rank 0), disagreeing with the 2D Delaunay
  // path; the rank only drives the early stop, so this was invisible in barcodes.
  std::vector<double> coords = {0.0, 0.0, 3.0, 4.0, 5.0, 0.0};
  Gudhi::reduced_rips::detail::Cloud pm{coords.data(), 2, 3};
  Gudhi::reduced_rips::Euclidean_kd_tree<> kd(pm, false);
  BOOST_CHECK_EQUAL(
      (Gudhi::reduced_rips::rng_cycle_rank_delaunay<double>(pm, kd, Gudhi::reduced_rips::delaunay_edges_2d)), 1u);
  BOOST_CHECK_EQUAL((Gudhi::reduced_rips::rng_cycle_rank_general<double>(pm, kd)), 1u);

  const double d12 = std::sqrt(20.0);
  Gudhi::reduced_rips::Matrix_geometry<double> geom({0.0, 5.0, 5.0, 5.0, 0.0, d12, 5.0, d12, 0.0}, 3);
  BOOST_CHECK_EQUAL(Gudhi::reduced_rips::rng_cycle_rank_matrix(geom), 1u);
}

BOOST_AUTO_TEST_CASE(integer_grid_ties_match_full_rips) {
  // An integer grid is saturated with exactly tied distances, including across the kd-tree's k-th-neighbor
  // boundary. The regression here is CGAL's arbitrary tie order (and tie-dependent k-subset) leaking into
  // the engine's heap frontier, which resumes positionally in a refreshed neighbor list: a mismatched
  // prefix skips one 1-simplex and doubles another, corrupting the filtration. Every unit square's loop is
  // a genuine bar (born 1, filled by the diagonal at sqrt 2), so both search strategies must reproduce the
  // full Vietoris-Rips ground truth exactly.
  Cloud grid;
  for (int x = 0; x < 5; ++x)
    for (int y = 0; y < 5; ++y) grid.push_back({double(x), double(y)});
  Bars truth = full_rips_h1(grid);
  for (auto strategy : {Reduced_rips::Search::kd_tree, Reduced_rips::Search::brute_force}) {
    Bars reduced = Reduced_rips::from_points(grid, 0, strategy).persistence();
    BOOST_CHECK(bars_close(truth, reduced));
  }
}

// ---- Degenerate and near-degenerate inputs -------------------------------------------------------------

BOOST_AUTO_TEST_CASE(coincident_points_do_not_break_the_barcode) {
  // A zero-distance pair of distinct points (here a duplicate of point 0) once crashed the point-cloud path:
  // the degenerate r == 0 lens query emitted a self-edge 2-simplex. Coincident points carry no H1, so the
  // barcode must match the duplicate-free circle, in both the 2D Delaunay and the >=4D general RNG paths.
  for (unsigned dim : {2u, 5u}) {
    BOOST_TEST_CONTEXT("dim=" << dim) {
      auto clean = circle(40, dim);
      auto dup = clean;
      dup.push_back(clean[0]);  // distance 0 between two distinct indices

      auto from_clean = Reduced_rips::from_points(clean);
      auto from_dup = Reduced_rips::from_points(dup);
      BOOST_REQUIRE(!from_clean.persistence().empty());  // the circle has H1; equal-but-empty is a failure
      BOOST_CHECK(bars_close(from_clean.persistence(), from_dup.persistence()));

      // The distance-matrix backend must agree on the same coincident-point input.
      Reduced_rips from_matrix = Reduced_rips::from_distance_matrix(full_distance_matrix(dup));
      BOOST_CHECK(bars_close(from_dup.persistence(), from_matrix.persistence()));
    }
  }
}

BOOST_AUTO_TEST_CASE(near_tie_distances_do_not_break_the_barcode) {
  // Regression for a floating-point hazard. A certified 2-simplex (a,b,c) records the ids of its three
  // edges, and the reduction relies on the two edges shorter than the diameter edge (a,b) having already been
  // popped and id'd. The index tie-break in the lune test enforces this for exact distances, but when two
  // edges are a near-tie in length, recomputing one of them with a different rounding (e.g. FMA contraction
  // under an optimized -march=native build) could admit a third vertex whose edge had not yet been id'd,
  // throwing std::out_of_range from the boundary lookup. A dense 3D torus has many such near-ties between
  // adjacent samples and reproduces the crash on an optimized build; the reduction must instead complete and
  // recover the torus's rank-2 H1 (two prominent loops). NOTE: this only exercises the hazard under
  // optimization with FMA contraction; an -O0 build cannot contract and so cannot reproduce it.
  Bars bc = Reduced_rips::from_points(torus(60, 40)).persistence();  // must not throw
  BOOST_CHECK_EQUAL(count_prominent(bc, 0.5), 2U);
  for (const auto& bar : bc) {
    BOOST_CHECK_GE(bar[0], 0.0);
    BOOST_CHECK_GT(bar[1], bar[0]);
    BOOST_CHECK(std::isfinite(bar[1]));
  }
}

// ---- Caching, diagnostics, integration, and empty/invalid inputs ---------------------------------------

BOOST_AUTO_TEST_CASE(persistence_is_cached_and_deterministic) {
  auto ph1 = Reduced_rips::from_points(circle(40, 2));
  auto first = ph1.persistence();  // copy
  const auto& second = ph1.persistence();
  BOOST_CHECK(first == second);

  // A fresh instance on identical input must give an identical barcode (no global RNG state leakage).
  auto other = Reduced_rips::from_points(circle(40, 2));
  BOOST_CHECK(other.persistence() == first);
}

BOOST_AUTO_TEST_CASE(diagnostic_counters_are_consistent) {
  auto ph1 = Reduced_rips::from_points(circle(40, 2));
  const auto& bc = ph1.persistence();
  // Recorded persistent pairs are exactly the (non-trivial) bars returned.
  BOOST_CHECK_EQUAL(ph1.num_persistence_pairs(), bc.size());
  // Every persistent pair needs a reduced 2-simplex column, and the loop produces at least one bar.
  BOOST_CHECK_GE(ph1.num_two_simplices(), ph1.num_persistence_pairs());
  BOOST_CHECK_GE(ph1.num_one_simplices(), 1u);
  BOOST_CHECK_GE(bc.size(), 1u);
}

BOOST_AUTO_TEST_CASE(off_file_point_cloud_integration) {
  // End-to-end on a real OFF file (a 3D torus of 300 points): exercises file I/O, the 3D Delaunay path at a
  // larger scale, and basic barcode well-formedness.
  Gudhi::Points_off_reader<std::vector<double>> off_reader("tore3D_300.off");
  BOOST_REQUIRE(off_reader.is_valid());
  auto ph1 = Reduced_rips::from_points(off_reader.get_point_cloud());
  const auto& bc = ph1.persistence();
  BOOST_REQUIRE(!bc.empty());
  for (const auto& bar : bc) {
    BOOST_CHECK_GE(bar[0], 0.0);
    BOOST_CHECK_GT(bar[1], bar[0]);
    BOOST_CHECK(std::isfinite(bar[1]));
  }
}

// Fewer than two points is not an error: the barcode is simply empty (there are no edges, hence no H1).
BOOST_AUTO_TEST_CASE(too_few_points_give_an_empty_barcode) {
  std::vector<std::vector<double>> no_points;
  BOOST_CHECK(Reduced_rips::from_points(no_points).persistence().empty());

  std::vector<std::vector<double>> one_point = {{0.0, 0.0}};
  auto from_one_point = Reduced_rips::from_points(one_point);
  BOOST_CHECK(from_one_point.persistence().empty());
  BOOST_CHECK_EQUAL(from_one_point.dimension(), 2U);

  // Zero-dimensional points carry no edges either: empty barcode, no error.
  std::vector<std::vector<double>> zero_dim = {{}, {}};
  BOOST_CHECK(Reduced_rips::from_points(zero_dim).persistence().empty());

  std::vector<std::vector<double>> empty_matrix;
  BOOST_CHECK(Reduced_rips::from_distance_matrix(empty_matrix).persistence().empty());

  std::vector<std::vector<double>> one_row = {{}};
  BOOST_CHECK(Reduced_rips::from_distance_matrix(one_row).persistence().empty());
}

BOOST_AUTO_TEST_CASE(degenerate_geometries_give_empty_barcodes) {
  // Collinear 2D points: the Delaunay triangulation is one-dimensional (edges but no faces), the RNG is the
  // path graph, and H1 is empty. Exercises the finite-edge extraction on a degenerate triangulation, where a
  // face-based extraction would find nothing and inflate the early-stop target.
  Cloud line;
  for (int i = 0; i < 20; ++i) line.push_back({double(i), 0.0});
  BOOST_CHECK(Reduced_rips::from_points(line).persistence().empty());

  // All points coincident: the Delaunay triangulation merges them into one vertex with no edges (cycle rank 0).
  Cloud same(10, {1.0, 2.0});
  BOOST_CHECK(Reduced_rips::from_points(same).persistence().empty());

  // Two points and a (scalene) triangle never carry H1.
  BOOST_CHECK(Reduced_rips::from_points(Cloud{{0.0, 0.0}, {1.0, 0.0}}).persistence().empty());
  BOOST_CHECK(Reduced_rips::from_points(Cloud{{0.0, 0.0}, {3.0, 0.0}, {3.0, 4.0}}).persistence().empty());
}

BOOST_AUTO_TEST_CASE(neighbor_budget_larger_than_n_is_clamped) {
  // A seed budget far beyond the point count must clamp to the available neighbors, not fail or change bars.
  auto pts = circle(12, 2);
  Bars budget_default = Reduced_rips::from_points(pts).persistence();
  Bars oversized = Reduced_rips::from_points(pts, 100).persistence();
  BOOST_REQUIRE(!budget_default.empty());
  BOOST_CHECK(bars_close(budget_default, oversized));
}

// The dimension-consistency check uses GUDHI_CHECK, which throws only in debug mode (GUDHI_DEBUG), so this is
// compiled only then. A genuine dimension mismatch among points is still rejected.
#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE(rejects_ragged_input) {
  std::vector<std::vector<double>> ragged = {{0.0, 0.0}, {1.0}};
  BOOST_CHECK_THROW(Reduced_rips::from_points(ragged), std::invalid_argument);
}
#endif
