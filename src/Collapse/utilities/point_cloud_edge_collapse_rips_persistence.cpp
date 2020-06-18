/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam, Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Flag_complex_edge_collapser.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/graph_simplicial_complex.h>

#include <boost/program_options.hpp>

#include<utility>  // for std::pair
#include<vector>

// Types definition

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Point = std::vector<Filtration_value>;
using Vector_of_points = std::vector<Point>;

using Flag_complex_edge_collapser = Gudhi::collapse::Flag_complex_edge_collapser<Vertex_handle, Filtration_value>;
using Proximity_graph = Gudhi::Proximity_graph<Flag_complex_edge_collapser>;

using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence);

int main(int argc, char* argv[]) {
  std::string off_file_points;
  std::string filediag;
  double threshold;
  int dim_max;
  int p;
  double min_persistence;

  program_options(argc, argv, off_file_points, filediag, threshold, dim_max, p, min_persistence);

  std::cout << "The current input values to run the program is: " << std::endl;
  std::cout << "min_persistence, threshold, max_complex_dimension, off_file_points, filediag"
            << std::endl;
  std::cout << min_persistence << ", " << threshold << ", " << dim_max
            << ", " << off_file_points << ", " << filediag << std::endl;

  Gudhi::Points_off_reader<Point> off_reader(off_file_points);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << off_file_points << "\n";
    exit(-1);  // ----- >>
  }

  Vector_of_points point_vector = off_reader.get_point_cloud();
  if (point_vector.size() <= 0) {
    std::cerr << "Empty point cloud." << std::endl;
    exit(-1);  // ----- >>
  }

  std::cout << "Successfully read " << point_vector.size() << " point_vector.\n";
  std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

  Proximity_graph proximity_graph = Gudhi::compute_proximity_graph<Simplex_tree>(point_vector,
                                                                                 threshold,
                                                                                 Gudhi::Euclidean_distance());

  if (num_edges(proximity_graph) <= 0) {
    std::cerr << "Total number of egdes are zero." << std::endl;
    exit(-1);
  }

  Flag_complex_edge_collapser edge_collapser(proximity_graph);

  Simplex_tree stree;
  for (Vertex_handle vertex = 0; static_cast<std::size_t>(vertex) < point_vector.size(); vertex++) {
    // insert the vertex with a 0. filtration value just like a Rips
    stree.insert_simplex({vertex}, 0.);
  }
  edge_collapser.process_edges(
    [&stree](Vertex_handle u, Vertex_handle v, Filtration_value filtration) {
        // insert the edge
        stree.insert_simplex({u, v}, filtration);
      });

  stree.expansion(dim_max);
  
  std::cout << "The complex contains " << stree.num_simplices() << " simplices  after collapse. \n";
  std::cout << "   and has dimension " << stree.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  stree.initialize_filtration();
  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(stree);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(p);

  pcoh.compute_persistent_cohomology(min_persistence);
  if (filediag.empty()) {
    pcoh.output_diagram();
  } else {
    std::ofstream out(filediag);
    pcoh.output_diagram(out);
    out.close();
  }

  return 0;
}

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
                       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")(
      "output-file,o", po::value<std::string>(&filediag)->default_value(std::string()),
      "Name of file in which the persistence diagram is written. Default print in std::cout")(
      "max-edge-length,r",
      po::value<Filtration_value>(&threshold)->default_value(std::numeric_limits<Filtration_value>::infinity()),
      "Maximal length of an edge for the Rips complex construction.")(
      "cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
      "Maximal dimension of the Rips complex we want to compute.")(
      "field-charac,p", po::value<int>(&p)->default_value(11),
      "Characteristic p of the coefficient field Z/pZ for computing homology.")(
      "min-persistence,m", po::value<Filtration_value>(&min_persistence),
      "Minimal lifetime of homology feature to be recorded. Default is 0. Enter a negative value to see zero length "
      "intervals");

  po::positional_options_description pos;
  pos.add("input-file", 1);

  po::options_description all;
  all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Compute the persistent homology with coefficient field Z/pZ \n";
    std::cout << "of a Rips complex, after edge collapse, defined on a set of input points.\n \n";
    std::cout << "The output diagram contains one bar per line, written with the convention: \n";
    std::cout << "   p   dim b d \n";
    std::cout << "where dim is the dimension of the homological feature,\n";
    std::cout << "b and d are respectively the birth and death of the feature and \n";
    std::cout << "p is the characteristic of the field Z/pZ used for homology coefficients." << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    exit(-1);
  }
}
