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

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_off_io.h>

#include <boost/program_options.hpp>

#include <string>
#include <vector>
#include <limits>  // infinity
#include <utility>  // for pair
#include <map>

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
using Proximity_graph = Gudhi::Proximity_graph<Graph>;

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
  Proximity_graph prox_graph = Gudhi::compute_proximity_graph<Graph>(off_reader.get_point_cloud(),
                                                                            threshold,
                                                                            Gudhi::Euclidean_distance());

  std::cout << "Vertices = " << boost::num_vertices(prox_graph) << std::endl;
  std::cout << "Edges = " << boost::num_edges(prox_graph) << std::endl;

  typename boost::graph_traits<Proximity_graph>::edge_iterator e_it, e_it_end;
  for (std::tie(e_it, e_it_end) = boost::edges(prox_graph); e_it != e_it_end;
       ++e_it) {
    auto u = source(*e_it, prox_graph);
    auto v = target(*e_it, prox_graph);
    if (v < u) std::swap(u, v);

    std::cout << "source = " << u << " - target = " << v <<  " - filtration = " << boost::get(Gudhi::edge_filtration_t(), prox_graph, *e_it) << std::endl;
  }

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
