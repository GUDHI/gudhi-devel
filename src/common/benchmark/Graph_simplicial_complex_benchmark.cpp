/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Clock.h>
#include <gudhi/Points_off_io.h>

#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for numeric limits
#include <fstream>
#include <cassert>


std::ofstream results_csv("results.csv");

template< typename Adjacency_list, typename ForwardPointRange, typename Distance >
Adjacency_list proximity_graph_computation(const ForwardPointRange& points, double threshold, Distance distance) {
  std::vector<std::pair< int, int >> edges;
  std::vector< double > edges_fil;
  std::map< int, double > vertices;

  int idx_u, idx_v;
  double fil;
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

  // Points are labeled from 0 to idx_u-1
  Adjacency_list skel_graph(edges.begin(), edges.end(), edges_fil.begin(), idx_u);

  auto vertex_prop = boost::get(Gudhi::vertex_filtration_t(), skel_graph);

  typename boost::graph_traits<Adjacency_list>::vertex_iterator vi, vi_end;
  for (std::tie(vi, vi_end) = boost::vertices(skel_graph);
       vi != vi_end; ++vi) {
    boost::put(vertex_prop, *vi, 0.);
  }

  return skel_graph;
}

template <typename Adjacency_list>
void benchmark_proximity_graph(const std::string& msg, const std::string& off_file_name) {
  Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file_name);
  assert(off_reader.is_valid());

  std::cout << "+ " << msg << std::endl;

  results_csv << "\"nb_points\";"
              << "\"nb_simplices\";"
              << "\"compute proximity graph(sec.)\";"
              << "\"complex_creation_time(sec.)\";"
              << "\"" << msg << "\";" << std::endl;

  Gudhi::Clock pg_compute_proximity_graph("    benchmark_proximity_graph - compute proximity graph");
  pg_compute_proximity_graph.begin();
  // benchmark begin
  Adjacency_list proximity_graph = proximity_graph_computation<Adjacency_list>(off_reader.get_point_cloud(),
                                                                               std::numeric_limits<double>::infinity(),
                                                                               Gudhi::Euclidean_distance());
  // benchmark end
  pg_compute_proximity_graph.end();
  std::cout << pg_compute_proximity_graph;

  Gudhi::Simplex_tree<> complex;
  Gudhi::Clock st_create_clock("    benchmark_proximity_graph - complex creation");
  st_create_clock.begin();
  // benchmark begin
  complex.insert_graph(proximity_graph);
  // benchmark end
  st_create_clock.end();
  std::cout << st_create_clock;

  results_csv << off_reader.get_point_cloud().size() << ";" << complex.num_simplices() << ";"
              << pg_compute_proximity_graph.num_seconds() << ";"
              << st_create_clock.num_seconds() << ";" << std::endl;

  std::cout << "    benchmark_proximity_graph - nb simplices = " << complex.num_simplices() << std::endl;
}

int main(int argc, char * const argv[]) {
  std::string off_file_name(argv[1]);

  // The fastest, the less memory used
  using vecSdirectedS = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                                              boost::property<Gudhi::vertex_filtration_t, double>,
                                              boost::property<Gudhi::edge_filtration_t, double>>;
  benchmark_proximity_graph<vecSdirectedS>("vecSdirectedS", off_file_name);

  using vecSundirectedS = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                                boost::property<Gudhi::vertex_filtration_t, double>,
                                                boost::property<Gudhi::edge_filtration_t, double>>;
  benchmark_proximity_graph<vecSundirectedS>("vecSundirectedS", off_file_name);

  using vecSbidirectionalS = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
                                                   boost::property<Gudhi::vertex_filtration_t, double>,
                                                   boost::property<Gudhi::edge_filtration_t, double>>;
  benchmark_proximity_graph<vecSbidirectionalS>("vecSbidirectionalS", off_file_name);

  using setSdirectedS = boost::adjacency_list<boost::setS, boost::vecS, boost::directedS,
                                              boost::property<Gudhi::vertex_filtration_t, double>,
                                              boost::property<Gudhi::edge_filtration_t, double>>;
  benchmark_proximity_graph<setSdirectedS>("setSdirectedS", off_file_name);

  using setSundirectedS = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                                boost::property<Gudhi::vertex_filtration_t, double>,
                                                boost::property<Gudhi::edge_filtration_t, double>>;
  benchmark_proximity_graph<setSundirectedS>("setSundirectedS", off_file_name);

  using setSbidirectionalS = boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS,
                                                   boost::property<Gudhi::vertex_filtration_t, double>,
                                                   boost::property<Gudhi::edge_filtration_t, double>>;
  benchmark_proximity_graph<setSbidirectionalS>("setSbidirectionalS", off_file_name);

  return 0;
}
