/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett, Musashi Koyama
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Delaunay_edges.h
 * @author Thomas Burnett, Musashi Koyama
 * @brief Edge extraction from the CGAL 2D/3D Delaunay triangulations, feeding the RNG cycle rank. This is the
 * only header of the module that uses the (GPL) CGAL Triangulation packages; it is kept separate so the rest
 * of the RNG machinery, and with it the whole distance-matrix path, stays free of CGAL includes.
 */

#ifndef REDUCED_RIPS_DELAUNAY_EDGES_H_
#define REDUCED_RIPS_DELAUNAY_EDGES_H_

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <gudhi/Reduced_rips/Helpers.h>

namespace Gudhi {

namespace reduced_rips {

// Deduplicated edges of the d-dimensional Delaunay triangulation, as pairs (first < second). Both dimensions
// read them from the finite edge iterator, which visits each edge exactly. Both edge sets are an Urquhart superset of
// the RNG.
inline std::vector<std::pair<std::size_t, std::size_t>> delaunay_edges_2d(const detail::Cloud& pm) {
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Vb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K>;
  using Tds = CGAL::Triangulation_data_structure_2<Vb>;
  using Delaunay = CGAL::Delaunay_triangulation_2<K, Tds>;
  using Point = K::Point_2;
  std::vector<std::pair<Point, std::size_t>> pts;
  pts.reserve(pm.size());
  for (std::size_t k = 0; k < pm.size(); ++k) pts.emplace_back(Point(pm[k][0], pm[k][1]), k);
  Delaunay t(pts.begin(), pts.end());
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  for (auto it = t.finite_edges_begin(); it != t.finite_edges_end(); ++it) {
    // A 2D edge is (face, i): the side of `face` opposite its vertex i, so its endpoints are the other two.
    std::size_t a = it->first->vertex(Delaunay::cw(it->second))->info();
    std::size_t b = it->first->vertex(Delaunay::ccw(it->second))->info();
    edges.emplace_back(std::minmax(a, b));
  }
  return edges;
}

inline std::vector<std::pair<std::size_t, std::size_t>> delaunay_edges_3d(const detail::Cloud& pm) {
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Vb = CGAL::Triangulation_vertex_base_with_info_3<std::size_t, K>;
  using Cb = CGAL::Delaunay_triangulation_cell_base_3<K>;
  using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
  using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>;
  using Point = Delaunay::Point;
  std::vector<std::pair<Point, std::size_t>> pts;
  pts.reserve(pm.size());
  for (std::size_t k = 0; k < pm.size(); ++k) pts.emplace_back(Point(pm[k][0], pm[k][1], pm[k][2]), k);
  Delaunay t(pts.begin(), pts.end());
  // The finite edge iterator visits each Delaunay edge exactly once, so we can emit
  // edges directly instead of deduplicating facet edges through a hash set.
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  edges.reserve(t.number_of_finite_edges());
  for (auto it = t.finite_edges_begin(); it != t.finite_edges_end(); ++it) {
    std::size_t a = it->first->vertex(it->second)->info();
    std::size_t b = it->first->vertex(it->third)->info();
    edges.emplace_back(std::minmax(a, b));
  }
  return edges;
}

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_DELAUNAY_EDGES_H_
