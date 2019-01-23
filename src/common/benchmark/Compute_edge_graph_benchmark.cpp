/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2019 Inria
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

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
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

struct Simplicial_complex {
  using Filtration_value = double;
  using Vertex_handle = int;
};

template<template<class> class Edge_graph, class SimplicialComplex>
void benchmark_compute_edge_graph(const std::string& msg, const std::string& off_file_name) {
  Gudhi::Points_off_reader<std::vector<typename SimplicialComplex::Filtration_value>> off_reader(off_file_name);
  assert(off_reader.is_valid());

  std::cout << "+ " << msg << std::endl;

  results_csv << "\"nb_points\";"
              << "\"compute edge graph(sec.)\";"
              << "\"" << msg << "\";" << std::endl;

  Gudhi::Clock pg_compute_edge_graph("    benchmark_compute_edge_graph - compute proximity graph");
  pg_compute_edge_graph.begin();
  // benchmark begin
  Edge_graph<SimplicialComplex> edge_graph =
      Gudhi::compute_edge_graph<Edge_graph, SimplicialComplex>(
          off_reader.get_point_cloud(),
          std::numeric_limits<typename SimplicialComplex::Filtration_value>::infinity(),
          Gudhi::Euclidean_distance());
  // benchmark end
  pg_compute_edge_graph.end();
  std::cout << pg_compute_edge_graph;

  results_csv << off_reader.get_point_cloud().size() << ";"
              << pg_compute_edge_graph.num_seconds() << ";" << std::endl;
}

template <typename SimplicialComplex>
using vecSdirectedS = typename boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS
    , boost::property < Gudhi::vertex_filtration_t, typename SimplicialComplex::Filtration_value >
    , boost::property < Gudhi::edge_filtration_t, typename SimplicialComplex::Filtration_value >>;

template <typename SimplicialComplex>
using vecSundirectedS = typename boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS
    , boost::property < Gudhi::vertex_filtration_t, typename SimplicialComplex::Filtration_value >
    , boost::property < Gudhi::edge_filtration_t, typename SimplicialComplex::Filtration_value >>;

template <typename SimplicialComplex>
using vecSbidirectionalS = typename boost::adjacency_list < boost::vecS, boost::vecS, boost::bidirectionalS
    , boost::property < Gudhi::vertex_filtration_t, typename SimplicialComplex::Filtration_value >
    , boost::property < Gudhi::edge_filtration_t, typename SimplicialComplex::Filtration_value >>;

template <typename SimplicialComplex>
using setSdirectedS = typename boost::adjacency_list < boost::setS, boost::vecS, boost::directedS
    , boost::property < Gudhi::vertex_filtration_t, typename SimplicialComplex::Filtration_value >
    , boost::property < Gudhi::edge_filtration_t, typename SimplicialComplex::Filtration_value >>;

template <typename SimplicialComplex>
using setSundirectedS = typename boost::adjacency_list < boost::setS, boost::vecS, boost::undirectedS
    , boost::property < Gudhi::vertex_filtration_t, typename SimplicialComplex::Filtration_value >
    , boost::property < Gudhi::edge_filtration_t, typename SimplicialComplex::Filtration_value >>;

template <typename SimplicialComplex>
using setSbidirectionalS = typename boost::adjacency_list < boost::setS, boost::vecS, boost::bidirectionalS
    , boost::property < Gudhi::vertex_filtration_t, typename SimplicialComplex::Filtration_value >
    , boost::property < Gudhi::edge_filtration_t, typename SimplicialComplex::Filtration_value >>;

int main(int argc, char * const argv[]) {
  std::string off_file_name(argv[1]);

  // The fastest, the less memory used - Not really efficient when sub-filtering is required (i.e. Strong collapse)
  benchmark_compute_edge_graph<vecSdirectedS, Simplicial_complex>("vecSdirectedS", off_file_name);

  benchmark_compute_edge_graph<vecSundirectedS, Simplicial_complex>("vecSundirectedS", off_file_name);

  benchmark_compute_edge_graph<vecSbidirectionalS, Simplicial_complex>("vecSbidirectionalS", off_file_name);

  /*benchmark_compute_edge_graph<setSdirectedS, Simplicial_complex>("setSdirectedS", off_file_name);

  benchmark_compute_edge_graph<setSundirectedS, Simplicial_complex>("setSundirectedS", off_file_name);

  benchmark_compute_edge_graph<setSbidirectionalS, Simplicial_complex>("setSbidirectionalS", off_file_name);
*/
  benchmark_compute_edge_graph<Gudhi::Filtered_edges_vector, Simplicial_complex>("Filtered_edges_vector", off_file_name);

  return 0;
}
