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

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using Graph_t = boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS
, boost::property < vertex_filtration_t, Filtration_value >
, boost::property < edge_filtration_t, Filtration_value >
>;
using Edge_t = std::pair< Vertex_handle, Vertex_handle >;

template< typename InputPointRange, typename Distance >
Graph_t compute_proximity_graph(InputPointRange &points, Filtration_value threshold, Distance distance);

using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;
using Point = std::vector<double>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;

void program_options(int argc, char * argv[]
                     , std::string & off_file_points
                     , std::string & filediag
                     , Filtration_value & threshold
                     , int & dim_max
                     , int & p
                     , Filtration_value & min_persistence);

int main(int argc, char * argv[]) {
  std::string off_file_points;
  std::string filediag;
  Filtration_value threshold;
  int dim_max;
  int p;
  Filtration_value min_persistence;

  program_options(argc, argv, off_file_points, filediag, threshold, dim_max, p, min_persistence);

  // Extract the points from the file filepoints
  Points_off_reader off_reader(off_file_points);

  // Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph(off_reader.get_point_cloud(), threshold
                                               , Euclidean_distance());

  // Construct the Rips complex in a Simplex Tree
  Simplex_tree st;
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion(dim_max);

  std::cout << "The complex contains " << st.num_simplices() << " simplices \n";
  std::cout << "   and has dimension " << st.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(st);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(p);

  pcoh.compute_persistent_cohomology(min_persistence);

  // Output the diagram in filediag
  if (filediag.empty()) {
    pcoh.output_diagram();
  } else {
    std::ofstream out(filediag);
    pcoh.output_diagram(out);
    out.close();
  }

  return 0;
}

void program_options(int argc, char * argv[]
                     , std::string & off_file_points
                     , std::string & filediag
                     , Filtration_value & threshold
                     , int & dim_max
                     , int & p
                     , Filtration_value & min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()
      ("input-file", po::value<std::string>(&off_file_points),
       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()
      ("help,h", "produce help message")
      ("output-file,o", po::value<std::string>(&filediag)->default_value(std::string()),
       "Name of file in which the persistence diagram is written. Default print in std::cout")
      ("max-edge-length,r",
       po::value<Filtration_value>(&threshold)->default_value(std::numeric_limits<Filtration_value>::infinity()),
       "Maximal length of an edge for the Rips complex construction.")
      ("cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
       "Maximal dimension of the Rips complex we want to compute.")
      ("field-charac,p", po::value<int>(&p)->default_value(11),
       "Characteristic p of the coefficient field Z/pZ for computing homology.")
      ("min-persistence,m", po::value<Filtration_value>(&min_persistence),
       "Minimal lifetime of homology feature to be recorded. Default is 0. Enter a negative value to see zero length intervals");

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
    std::cout << "Compute the persistent homology with coefficient field Z/pZ \n";
    std::cout << "of a Rips complex defined on a set of input points.\n \n";
    std::cout << "The output diagram contains one bar per line, written with the convention: \n";
    std::cout << "   p   dim b d \n";
    std::cout << "where dim is the dimension of the homological feature,\n";
    std::cout << "b and d are respectively the birth and death of the feature and \n";
    std::cout << "p is the characteristic of the field Z/pZ used for homology coefficients." << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}

/** Output the proximity graph of the points.
 *
 * If points contains n elements, the proximity graph is the graph
 * with n vertices, and an edge [u,v] iff the distance function between
 * points u and v is smaller than threshold.
 *
 * The type PointCloud furnishes .begin() and .end() methods, that return
 * iterators with value_type Point.
 */
template< typename InputPointRange, typename Distance >
Graph_t compute_proximity_graph(InputPointRange &points, Filtration_value threshold, Distance distance) {
  std::vector< Edge_t > edges;
  std::vector< Filtration_value > edges_fil;

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
