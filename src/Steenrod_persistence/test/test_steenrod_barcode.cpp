/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Anibal M. Medina-Mardones
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <algorithm>
#include <set>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "steenrod_persistence"
#include <boost/test/unit_test.hpp>

#include <gudhi/Steenrod_barcode.h>

using namespace Gudhi::steenrod_persistence;

BOOST_AUTO_TEST_CASE(symm_diff_basic) {
  BOOST_CHECK((symm_diff({1, 3}, {2, 4}) == Column{1, 2, 3, 4}));
  BOOST_CHECK((symm_diff({1, 2, 3}, {2, 3, 4}) == Column{1, 4}));
  BOOST_CHECK(symm_diff({1, 2, 3}, {1, 2, 3}).empty());
  BOOST_CHECK((symm_diff({}, {1, 2}) == Column{1, 2}));
}

BOOST_AUTO_TEST_CASE(pivot_basic) {
  BOOST_CHECK_EQUAL(pivot(Column{}), -1);
  BOOST_CHECK_EQUAL(pivot(Column{3, 7, 9}), 3);
}

// Build a Filtration_by_dim from a list of top simplices — closes the complex
// by taking all subsets and assigns filtration values in insertion order.
static Filtration_by_dim make_filtration_from_top(const std::vector<Simplex>& top) {
  // 1. Close the complex: collect every subset of every top simplex.
  std::set<Simplex> all_simplices;
  for (const Simplex& t : top) {
    const int n = static_cast<int>(t.size());
    for (int mask = 1; mask < (1 << n); ++mask) {
      Simplex s;
      for (int i = 0; i < n; ++i) {
        if (mask & (1 << i)) s.push_back(t[i]);
      }
      std::sort(s.begin(), s.end());
      all_simplices.insert(s);
    }
  }

  // 2. Sort simplices by (dim, lexicographic) — dim-first ordering is the
  //    simplest filtration that closes downward correctly.
  std::vector<Simplex> sorted(all_simplices.begin(), all_simplices.end());
  std::sort(sorted.begin(), sorted.end(), [](const Simplex& a, const Simplex& b) {
    if (a.size() != b.size()) return a.size() < b.size();
    return a < b;
  });

  // 3. Bucket by dimension with absolute indices in insertion order.
  Filtration_by_dim fbd;
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    const int dim = static_cast<int>(sorted[i].size()) - 1;
    while (static_cast<int>(fbd.size()) <= dim) fbd.emplace_back();
    fbd[dim].idxs.push_back(static_cast<Index>(i));
    fbd[dim].tups.push_back(sorted[i]);
  }
  return fbd;
}

// Minimal triangulation of the real projective plane RP^2 (6 vertices,
// 15 edges, 10 triangles). Its Sq^1 barcode has exactly one essential bar in
// homological dimension 2.
BOOST_AUTO_TEST_CASE(rp2_steenrod_sq1_essential) {
  const std::vector<Simplex> top = {
      {1, 2, 4}, {2, 3, 4}, {1, 3, 5}, {2, 3, 5}, {1, 4, 5},
      {1, 2, 6}, {1, 3, 6}, {3, 4, 6}, {2, 5, 6}, {4, 5, 6}};
  const auto fbd = make_filtration_from_top(top);

  const auto result = barcodes(/*k=*/1, fbd);

  // Ordinary barcode sanity: RP^2 over F_2 is acyclic-free —
  // H^0 = F_2, H^1 = F_2, H^2 = F_2 — so three essential bars total.
  int essential = 0;
  for (const auto& bars_dim : result.ordinary) {
    for (const auto& bar : bars_dim) {
      if (bar[0] == -1) ++essential;
    }
  }
  BOOST_CHECK_EQUAL(essential, 3);

  // Sq^1 essential bar in dim 2: detects the non-orientability of RP^2.
  BOOST_REQUIRE_GE(static_cast<int>(result.steenrod.size()), 3);
  int sq1_essential_dim2 = 0;
  for (const auto& bar : result.steenrod[2]) {
    if (bar[0] == -1) ++sq1_essential_dim2;
  }
  BOOST_CHECK_EQUAL(sq1_essential_dim2, 1);
}
