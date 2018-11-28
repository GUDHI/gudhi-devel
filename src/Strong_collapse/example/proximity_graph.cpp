/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
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

#include <gudhi/distance_functions.h>
#include <gudhi/Points_off_io.h>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_for.h>
#endif

#include <boost/program_options.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <utility>  // for pair<>
#include <vector>
#include <map>
#include <tuple>  // for std::tie

#include <string>
#include <vector>
#include <limits>  // infinity
#include <utility>  // for pair

/* Edge tag for Boost PropertyGraph. */
struct edge_filtration_t {
  typedef boost::edge_property_tag kind;
};

/** \brief Proximity_graph contains the edges with their filtration values and vertices in order to store the result
 * of `compute_proximity_graph` function.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Filtration_value` type definition.
 *
 */
template <typename SimplicialComplexForProximityGraph>
using Proximity_graph = typename boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property,  // no edge property
    boost::property < edge_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value >>;

// setS returns pointer-like vertices (don't know why):
// source = 0x55ef4a232310 - target = 0x55ef4a232670 - filtration = 6.08276
// instead of :
// source = 0 - target = 1 - filtration = 6.08276
// template <typename SimplicialComplexForProximityGraph>
// using Proximity_graph = typename boost::adjacency_list < boost::setS, boost::setS, boost::undirectedS
//     , boost::property < vertex_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value >
//     , boost::property < edge_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value >>;

template< typename SimplicialComplexForProximityGraph
    , typename ForwardPointRange
    , typename Distance >
Proximity_graph<SimplicialComplexForProximityGraph> compute_proximity_graph(
    const ForwardPointRange& points,
    typename SimplicialComplexForProximityGraph::Filtration_value threshold,
    Distance distance) {
  using Vertex_handle = typename SimplicialComplexForProximityGraph::Vertex_handle;
  using Filtration_value = typename SimplicialComplexForProximityGraph::Filtration_value;

  std::vector<std::pair< Vertex_handle, Vertex_handle >> edges;
  std::vector< Filtration_value > edges_fil;

  Vertex_handle idx_u, idx_v;
  Filtration_value fil;
  idx_u = 0;
  for (auto it_u = points.begin(); it_u != points.end(); ++it_u) {
    idx_v = idx_u + 1;
    for (auto it_v = it_u + 1; it_v != points.end(); ++it_v, ++idx_v) {
      fil = distance(*it_u, *it_v);
      if (fil <= threshold) {
        auto edges_iter = edges.begin();
        auto edges_fil_iter = edges_fil.begin();
        while (edges_iter < edges.end() && edges_fil_iter != edges_fil.end() && fil > *edges_fil_iter) {
          ++edges_iter;
          ++edges_fil_iter;
        }

        edges.insert(edges_iter, std::make_pair(idx_u, idx_v));
        edges_fil.insert(edges_fil_iter, fil);
      }
    }
    ++idx_u;
  }

  // Points are labeled from 0 to idx_u-1
  Proximity_graph<SimplicialComplexForProximityGraph> skel_graph(edges.begin(), edges.end(), edges_fil.begin(), idx_u);

  return skel_graph;
}

// ----------------------------------------------------------------------------
// rips_persistence_step_by_step is an example of each step that is required to
// build a Rips over a Simplex_tree. Please refer to rips_persistence to see
// how to do the same thing with the Rips_complex wrapper for less detailed
// steps.
// ----------------------------------------------------------------------------

// Types definition
struct Graph {
  using Vertex_handle = int;
  using Filtration_value = double;
};

using Vertex_handle = Graph::Vertex_handle;
using Filtration_value = Graph::Filtration_value;

using Point = std::vector<double>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;

void program_options(int argc, char * argv[]
    , std::string & off_file_points
    , Filtration_value & threshold);

int main(int argc, char * argv[]) {
  std::string off_file_points;
  Filtration_value threshold;

  program_options(argc, argv, off_file_points, threshold);

  // Extract the points from the file filepoints
  Points_off_reader off_reader(off_file_points);

  // Compute the proximity graph of the points
  Proximity_graph<Graph> prox_graph = compute_proximity_graph<Graph>(off_reader.get_point_cloud(),
                                                                            threshold,
                                                                            Gudhi::Euclidean_distance());

  std::cout << "Vertices = " << boost::num_vertices(prox_graph) << std::endl;
  std::cout << "Edges = " << boost::num_edges(prox_graph) << std::endl;

  std::vector<Filtration_value> alpha_values;
  for (std::vector<Filtration_value>::size_type index = 0; index < 10; index ++) {
    alpha_values.push_back((index + 1) * threshold/10.);  // Can be 0
  }

#if defined(GUDHI_USE_TBB)
  tbb::parallel_for<std::size_t>( static_cast<std::size_t>(0), alpha_values.size(), [&](std::size_t index){
    Filtration_value alpha = alpha_values[index];
    typename boost::graph_traits<Proximity_graph<Graph>>::edge_iterator e_it, e_it_end;
    for (std::tie(e_it, e_it_end) = boost::edges(prox_graph); e_it != e_it_end;
         ++e_it) {
      auto u = source(*e_it, prox_graph);
      auto v = target(*e_it, prox_graph);
      //if (v < u) std::swap(u, v);
      if (boost::get(edge_filtration_t(), prox_graph, *e_it) <= alpha) {
        std::cout << "alpha = " << alpha << " - source = " << u << " - target = " << v << " - filtration = "
                  << boost::get(edge_filtration_t(), prox_graph, *e_it) << std::endl;
      } else {
        break;
      }
    }
  });
#else
  for (Filtration_value alpha : alpha_values) {
    typename boost::graph_traits<Proximity_graph<Graph>>::edge_iterator e_it, e_it_end;
    for (std::tie(e_it, e_it_end) = boost::edges(prox_graph); e_it != e_it_end;
         ++e_it) {
      auto u = source(*e_it, prox_graph);
      auto v = target(*e_it, prox_graph);
      //if (v < u) std::swap(u, v);
      if (boost::get(edge_filtration_t(), prox_graph, *e_it) <= alpha) {
        std::cout << "alpha = " << alpha << "source = " << u << " - target = " << v << " - filtration = "
                  << boost::get(edge_filtration_t(), prox_graph, *e_it) << std::endl;
      } else {
        break;
      }
    }
  }
#endif

  return 0;
}

void program_options(int argc, char * argv[]
    , std::string & off_file_points
    , Filtration_value & threshold) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()
      ("input-file", po::value<std::string>(&off_file_points),
       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 5);
  visible.add_options()
      ("help,h", "produce help message")
      ("max-edge-length,r",
       po::value<Filtration_value>(&threshold)->default_value(std::numeric_limits<Filtration_value>::infinity()),
       "Maximal length of an edge for the Rips complex construction.");

  po::positional_options_description pos;
  pos.add("input-file", 1);

  po::options_description all;
  all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
      options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Compute the proximity graph on a set of input points." << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    exit(-1);
  }
}
