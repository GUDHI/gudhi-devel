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

#ifndef UTILS_EDGE_COLLAPSOR_H_
#define UTILS_EDGE_COLLAPSOR_H_

#include <list>
#include "utils/Edge_contractor.h"
#include "utils/UI_utils.h"

/**
 * Iteratively puts every vertex at the center of its neighbors
 */
template<typename SkBlComplex> class Edge_collapsor {
 private:
  SkBlComplex& complex_;
  unsigned num_collapses_;

 public:
  typedef typename SkBlComplex::Vertex_handle Vertex_handle;
  typedef typename SkBlComplex::Edge_handle Edge_handle;

  /**
   * @brief Collapse num_collapses edges. If num_collapses<0 then it collapses all possible edges.
   * Largest edges are collapsed first.
   *
   */
  Edge_collapsor(SkBlComplex& complex, unsigned num_collapses) :
      complex_(complex), num_collapses_(num_collapses) {
    std::list<Edge_handle> edges;
    edges.insert(edges.begin(), complex_.edge_range().begin(), complex_.edge_range().end());

    edges.sort(
               [&](Edge_handle e1, Edge_handle e2) {
                 return squared_edge_length(e1) < squared_edge_length(e2);
               });

    collapse_edges(edges);
  }

 private:
  void collapse_edges(std::list<Edge_handle>& edges) {
    while (!edges.empty() && num_collapses_--) {
      Edge_handle current_edge = edges.front();
      edges.pop_front();
      if (is_link_reducible(current_edge))
        complex_.remove_edge(current_edge);
    }
  }

  bool is_link_reducible(Edge_handle e) {
    auto link = complex_.link(e);

    if (link.empty())
      return false;

    if (link.is_cone())
      return true;

    if (link.num_connected_components() > 1)
      return false;

    Edge_contractor<SkBlComplex> contractor(link, link.num_vertices() - 1);

    return (link.num_vertices() == 1);
  }

  double squared_edge_length(Edge_handle e) const {
    return squared_eucl_distance(complex_.point(complex_.first_vertex(e)), complex_.point(complex_.second_vertex(e)));
  }

  double squared_eucl_distance(const Point& p1, const Point& p2) const {
    return Geometry_trait::Squared_distance_d()(p1, p2);
  }
};

#endif  // UTILS_EDGE_COLLAPSOR_H_
