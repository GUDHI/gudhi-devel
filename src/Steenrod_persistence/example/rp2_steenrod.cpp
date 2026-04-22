/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Maximiliano Alvarez
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include <gudhi/Steenrod_barcode.h>

using namespace Gudhi::steenrod_persistence;

// Close the complex spanned by ``top`` and return a Filtration_by_dim whose
// filtration order is (dim, lexicographic).
static Filtration_by_dim make_filtration_from_top(const std::vector<Simplex>& top) {
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

  std::vector<Simplex> sorted(all_simplices.begin(), all_simplices.end());
  std::sort(sorted.begin(), sorted.end(), [](const Simplex& a, const Simplex& b) {
    if (a.size() != b.size()) return a.size() < b.size();
    return a < b;
  });

  Filtration_by_dim fbd;
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    const int dim = static_cast<int>(sorted[i].size()) - 1;
    while (static_cast<int>(fbd.size()) <= dim) fbd.emplace_back();
    fbd[dim].idxs.push_back(static_cast<Index>(i));
    fbd[dim].tups.push_back(sorted[i]);
  }
  return fbd;
}

// Minimal triangulation of the real projective plane.
// 6 vertices, 15 edges, 10 triangles.
int main() {
  const std::vector<Simplex> top = {
      {1, 2, 4}, {2, 3, 4}, {1, 3, 5}, {2, 3, 5}, {1, 4, 5},
      {1, 2, 6}, {1, 3, 6}, {3, 4, 6}, {2, 5, 6}, {4, 5, 6}};
  const auto fbd = make_filtration_from_top(top);
  const auto result = barcodes(/*k=*/1, fbd);

  std::cout << "RP^2 — Sq^1 Steenrod barcode, bars per dimension:\n";
  for (std::size_t d = 0; d < result.steenrod.size(); ++d) {
    std::cout << "  dim " << d << " : " << result.steenrod[d].size() << " bar(s)\n";
  }
  return 0;
}
