/* This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#ifndef UTILS_CRITICAL_POINTS_H_
#define UTILS_CRITICAL_POINTS_H_

#include <deque>
#include <utility>  // for pair<>
#include <algorithm>  // for sort

#include "utils/Edge_contractor.h"

/**
 * Iteratively tries to anticollapse smallest edge non added so far.
 * If its link is contractible then no topological change and else possible topological change.
 *
 * todo do a sparsification with some parameter eps while growing
 */
template<typename SkBlComplex> class Critical_points {
 private:
  SkBlComplex filled_complex_;
  const SkBlComplex& input_complex_;
  double max_length_;
  std::ostream& stream_;

 public:
  typedef typename SkBlComplex::Vertex_handle Vertex_handle;
  typedef typename SkBlComplex::Edge_handle Edge_handle;
  typedef typename std::pair<Vertex_handle, Vertex_handle> Edge;

  /**
   * @brief check all pair of points with length smaller than max_length
   */
  Critical_points(const SkBlComplex& input_complex, std::ostream& stream, double max_length) :
      input_complex_(input_complex), max_length_(max_length), stream_(stream) {
    std::deque<Edge> edges;
    auto vertices = input_complex.vertex_range();
    for (auto p = vertices.begin(); p != vertices.end(); ++p) {
      filled_complex_.add_vertex(input_complex.point(*p));
      for (auto q = p; ++q != vertices.end(); /**/)
        if (squared_eucl_distance(input_complex.point(*p), input_complex.point(*q)) < max_length_ * max_length_)
          edges.emplace_back(*p, *q);
    }

    std::sort(edges.begin(), edges.end(),
              [&](Edge e1, Edge e2) {
                return squared_edge_length(e1) < squared_edge_length(e2);
              });

    anti_collapse_edges(edges);
  }

 private:
  double squared_eucl_distance(const Point& p1, const Point& p2) const {
    return Geometry_trait::Squared_distance_d()(p1, p2);
  }

  void anti_collapse_edges(const std::deque<Edge>& edges) {
    unsigned pos = 0;
    for (Edge e : edges) {
      std::cout << "edge " << pos++ << "/" << edges.size() << "\n";
      auto eh = filled_complex_.add_edge(e.first, e.second);
      int is_contractible(is_link_reducible(eh));

      switch (is_contractible) {
        case 0:
          stream_ << "alpha=" << std::sqrt(squared_edge_length(e)) << " topological change" << std::endl;
          break;
        case 2:
          stream_ << "alpha=" << std::sqrt(squared_edge_length(e)) << " maybe a topological change" << std::endl;
          break;
        default:
          break;
      }
    }
  }

  // 0 -> not
  // 1 -> yes
  // 2 -> maybe

  int is_link_reducible(Edge_handle e) {
    auto link = filled_complex_.link(e);

    if (link.empty())
      return 0;

    Edge_contractor<Complex> contractor(link, link.num_vertices() - 1);
    (void)contractor;

    if (link.num_connected_components() > 1)
      // one than more CC -> not contractible
      return 0;

    if (link.num_vertices() == 1)
      // reduced to one point -> contractible
      return 1;
    else
      // we dont know
      return 2;
  }

  double squared_edge_length(Edge_handle e) const {
    return squared_eucl_distance(input_complex_.point(input_complex_.first_vertex(e)),
                                 input_complex_.point(input_complex_.second_vertex(e)));
  }

  double squared_edge_length(Edge e) const {
    return squared_eucl_distance(input_complex_.point(e.first), input_complex_.point(e.second));
  }
};

#endif  // UTILS_CRITICAL_POINTS_H_
