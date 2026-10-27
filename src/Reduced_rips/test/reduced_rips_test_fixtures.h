/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

// Shared point-cloud generators, barcode helpers and the full-Vietoris-Rips ground-truth pipeline for the
// Reduced_rips unit tests. Each test file defines its own BOOST_TEST_MODULE and then includes this header.

#ifndef REDUCED_RIPS_TEST_FIXTURES_H_
#define REDUCED_RIPS_TEST_FIXTURES_H_

#include <gudhi/Reduced_rips.h>

#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/distance_functions.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <utility>
#include <vector>

using Reduced_rips = Gudhi::reduced_rips::Reduced_rips<>;
// A barcode as an array of {birth, death} pairs of doubles: the shape Reduced_rips<>::persistence() returns.
using Bars = std::vector<std::array<double, 2>>;
using Cloud = std::vector<std::vector<double>>;

// Standard full-complex degree-1 PH pipeline, used as ground truth (with Z/2Z coefficients, matching
// Reduced_rips). A double-precision Simplex_tree keeps the comparison tight.
using Stree = Gudhi::Simplex_tree<>;
using Filtration_value = Stree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Stree, Field_Zp>;

// ---- Point-cloud generators ----------------------------------------------------------------------------

inline constexpr double pi = 3.14159265358979323846;

// n points evenly sampled on the unit circle in dimension `dim` (padded with zeros). Degree-1 PH of such
// a sample has exactly one prominent class (the loop), born around the sampling spacing and dying near the
// circle diameter.
inline Cloud circle(unsigned n, std::size_t dim = 2) {
  Cloud pts;
  for (unsigned i = 0; i < n; ++i) {
    std::vector<double> p(dim, 0.0);
    double theta = 2.0 * pi * i / n;
    p[0] = std::cos(theta);
    p[1] = std::sin(theta);
    pts.push_back(std::move(p));
  }
  return pts;
}

// Two unit circles far apart on the x-axis: two independent H1 loops.
inline Cloud two_circles(unsigned n) {
  Cloud a = circle(n, 2);
  Cloud b = circle(n, 2);
  for (auto& p : b) p[0] += 10.0;
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

// A torus in R^3 (nu points around the tube times nv around it): H1 of rank 2. Small enough that the full
// Rips ground truth is cheap, large enough to exercise the 3D Delaunay relative-neighborhood-graph path.
inline Cloud torus(unsigned nu, unsigned nv, double R = 2.0, double r = 0.7) {
  Cloud pts;
  for (unsigned i = 0; i < nu; ++i) {
    double u = 2.0 * pi * i / nu;
    for (unsigned j = 0; j < nv; ++j) {
      double v = 2.0 * pi * j / nv;
      pts.push_back({(R + (r * std::cos(v))) * std::cos(u), (R + (r * std::cos(v))) * std::sin(u), r * std::sin(v)});
    }
  }
  return pts;
}

// Deterministic uniform random cloud in the unit cube of the given dimension.
inline Cloud random_cloud(unsigned n, std::size_t dim, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  Cloud pts(n, std::vector<double>(dim));
  for (auto& p : pts)
    for (auto& c : p) c = unit(gen);
  return pts;
}

// Full (symmetric, n x n) Euclidean distance matrix of a point cloud.
inline Cloud full_distance_matrix(const Cloud& pts) {
  std::size_t n = pts.size();
  Cloud m(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      double s = 0.0;
      for (std::size_t k = 0; k < pts[i].size(); ++k) {
        double d = pts[i][k] - pts[j][k];
        s += d * d;
      }
      m[i][j] = std::sqrt(s);
    }
  return m;
}

// ---- Barcode helpers -----------------------------------------------------------------------------------

// Finite bars of positive persistence, sorted. Reduced_rips already drops zero-length bars; applying the
// same filter to both sides makes the comparison robust to near-zero numerical noise.
inline Bars finite_positive(const Bars& b, double eps = 1e-9) {
  Bars out;
  for (const auto& x : b)
    if (std::isfinite(x[1]) && x[1] - x[0] > eps) out.push_back(x);
  std::sort(out.begin(), out.end());
  return out;
}

// True iff two barcodes match as multisets within tolerance. Distances squared then square-rooted along the
// matrix path can differ from the coordinate path by a few ULPs, hence the tolerance rather than equality.
inline bool bars_close(Bars a, Bars b, double tol = 1e-6) {
  a = finite_positive(a);
  b = finite_positive(b);
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (std::abs(a[i][0] - b[i][0]) > tol || std::abs(a[i][1] - b[i][1]) > tol) return false;
  return true;
}

// Number of bars more persistent than `min_persistence` (counts prominent topological features).
inline std::size_t count_prominent(const Bars& b, double min_persistence) {
  std::size_t c = 0;
  for (const auto& x : b)
    if (x[1] - x[0] > min_persistence) ++c;
  return c;
}

// Re-express a Filtration_value-typed barcode as the double-valued Bars the comparison helpers consume.
template <class FV>
Bars to_double_bars(const std::vector<std::array<FV, 2>>& b) {
  Bars out;
  out.reserve(b.size());
  for (const auto& x : b) out.push_back({static_cast<double>(x[0]), static_cast<double>(x[1])});
  return out;
}

// Ground-truth degree-1 barcode of the full (un-thresholded) Vietoris-Rips filtration.
inline Bars full_rips_h1(const Cloud& pts) {
  Rips_complex rips(pts, std::numeric_limits<double>::infinity(), Gudhi::Euclidean_distance());
  Stree st;
  rips.create_complex(st, 2);
  Persistent_cohomology pcoh(st);
  pcoh.init_coefficients(2);
  pcoh.compute_persistent_cohomology(0.0);
  Bars bars;
  for (const auto& bd : pcoh.intervals_in_dimension(1)) bars.push_back({bd.first, bd.second});
  return bars;
}

// Ground-truth degree-1 barcode of the full Vietoris-Rips filtration of a (full, symmetric) distance matrix.
// 2-simplices fill at their longest edge, matching the Reduced_rips convention; no triangle inequality needed.
inline Bars full_rips_h1_from_matrix(const std::vector<std::vector<double>>& full) {
  Stree st;
  for (std::size_t i = 0; i < full.size(); ++i)
    for (std::size_t j = i + 1; j < full.size(); ++j)
      st.insert_simplex_and_subfaces({static_cast<int>(i), static_cast<int>(j)}, full[i][j]);
  st.expansion(2);
  st.make_filtration_non_decreasing();
  Persistent_cohomology pcoh(st);
  pcoh.init_coefficients(2);
  pcoh.compute_persistent_cohomology(0.0);
  Bars truth;
  for (const auto& bd : pcoh.intervals_in_dimension(1)) truth.push_back({bd.first, bd.second});
  return truth;
}

#endif  // REDUCED_RIPS_TEST_FIXTURES_H_
