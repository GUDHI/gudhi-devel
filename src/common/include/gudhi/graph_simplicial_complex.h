/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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
 */

#ifndef GRAPH_SIMPLICIAL_COMPLEX_H_
#define GRAPH_SIMPLICIAL_COMPLEX_H_

#include <boost/graph/adjacency_list.hpp>

#include <utility>  // for pair<>
#include <vector>
#include <map>

/* Edge tag for Boost PropertyGraph. */
struct edge_filtration_t {
  typedef boost::edge_property_tag kind;
};

/* Vertex tag for Boost PropertyGraph. */
struct vertex_filtration_t {
  typedef boost::vertex_property_tag kind;
};

typedef int Vertex_handle;
typedef double Filtration_value;
typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS
, boost::property < vertex_filtration_t, Filtration_value >
, boost::property < edge_filtration_t, Filtration_value >
> Graph_t;
typedef std::pair< Vertex_handle, Vertex_handle > Edge_t;

/** \brief Output the proximity graph of the points.
 *
 * If points contains n elements, the proximity graph is the graph 
 * with n vertices, and an edge [u,v] iff the distance function between 
 * points u and v is smaller than threshold.
 *
 * The type PointCloud furnishes .begin() and .end() methods, that return
 * iterators with value_type Point.
 */
template< typename PointCloud
, typename Point >
Graph_t compute_proximity_graph(PointCloud &points
                                , Filtration_value threshold
                                , Filtration_value distance(Point p1, Point p2)) {
  std::vector< Edge_t > edges;
  std::vector< Filtration_value > edges_fil;
  std::map< Vertex_handle, Filtration_value > vertices;

  Vertex_handle idx_u, idx_v;
  Filtration_value fil;
  idx_u = 0;
  for (auto it_u = points.begin(); it_u != points.end(); ++it_u) {
    idx_v = idx_u + 1;
    for (auto it_v = it_u + 1; it_v != points.end(); ++it_v, ++idx_v) {
      fil = distance(*it_u, *it_v);
      if (fil <= threshold) {
        edges.emplace_back(idx_u, idx_v);
        edges_fil.push_back(fil);
      }
    }
    ++idx_u;
  }

  Graph_t skel_graph(edges.begin()
                     , edges.end()
                     , edges_fil.begin()
                     , idx_u);  // number of points labeled from 0 to idx_u-1

  auto vertex_prop = boost::get(vertex_filtration_t(), skel_graph);

  boost::graph_traits<Graph_t>::vertex_iterator vi, vi_end;
  for (std::tie(vi, vi_end) = boost::vertices(skel_graph);
       vi != vi_end; ++vi) {
    boost::put(vertex_prop, *vi, 0.);
  }

  return skel_graph;
}

#endif  // GRAPH_SIMPLICIAL_COMPLEX_H_
