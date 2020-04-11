/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam, Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Flag_complex_sparse_matrix.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Rips_edge_list.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Points_off_io.h>

#include <boost/program_options.hpp>

using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = Simplex_tree::Filtration_value;
using Vertex_handle = Simplex_tree::Vertex_handle;

using Flag_complex_sparse_matrix = Gudhi::collapse::Flag_complex_sparse_matrix<Vertex_handle, Filtration_value>;

using Rips_edge_list = Gudhi::rips_edge_list::Rips_edge_list<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
using Distance_matrix = std::vector<std::vector<Filtration_value>>;

void program_options(int argc, char* const argv[], double& min_persistence, double& end_thresold,
                     int& dimension, int& dim_max, std::string& csv_matrix_file, std::string& filediag) {
  namespace po = boost::program_options;
  po::options_description visible("Allowed options", 100);
      visible.add_options()
        ("help,h", "produce help message")
      	("min_persistence,m", po::value<double>(&min_persistence)->default_value(0.1),
         "Minimum persistence interval length")
        ("end_thresold,e", po::value<double>(&end_thresold)->default_value(1),
         "Final threshold for rips complex.")
        ("dimensions,D", po::value<int>(&dimension)->default_value(2),
         "Dimension of the manifold.")
        ("dim_max,k ", po::value<int>(&dim_max)->default_value(2),
         "Maximum allowed dimension of the Rips complex.")
        ("input_file_name,i", po::value<std::string>(&csv_matrix_file),
         "The input file.")
        ("filediag,o", po::value<std::string>(&filediag),
         "The output file.");

  po::options_description all;
  all.add(visible);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << std::endl;
    std::cout << "Computes rips complexes of different threshold values, to 'end_thresold' from a n random uniform "
                 "point_vector on a selected manifold, . \n";
    std::cout << "Strongly collapses all the rips complexes and output the results in out_file. \n";
    std::cout << "The experiments are repeted 'repete' num of times for each threshold value. \n";
    std::cout << "type -m for manifold options, 's' for uni sphere, 'b' for unit ball, 'f' for file. \n";
    std::cout << "type -i 'filename' for Input file option for exported point sample. \n";
    std::cout << std::endl << std::endl;
    std::cout << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}

void program_options(int argc, char* argv[], std::string& csv_matrix_file, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence);

int main(int argc, char* argv[]) {
  std::string csv_matrix_file;
  std::string filediag;
  Filtration_value threshold;
  int dim_max = 2;
  int p;
  Filtration_value min_persistence;

  program_options(argc, argv, csv_matrix_file, filediag, threshold, dim_max, p, min_persistence);

  Distance_matrix distances;
  Distance_matrix sparse_distances;

  distances = Gudhi::read_lower_triangular_matrix_from_csv_file<Filtration_value>(csv_matrix_file);
  std::cout << "Read the distance matrix succesfully, of size: " << distances.size() << std::endl;

  Flag_complex_sparse_matrix::Filtered_sorted_edge_list edge_t;

  Rips_edge_list Rips_edge_list_from_file(distances, threshold);
  Rips_edge_list_from_file.create_edges(edge_t);
  
  if (edge_t.size() <= 0) {
    std::cerr << "Total number of egdes are zero." << std::endl;
    exit(-1);
  }

  std::cout << "Total number of edges before collapse are: " << edge_t.size() << std::endl;

  // Now we will perform filtered edge collapse to sparsify the edge list edge_t.
  Flag_complex_sparse_matrix mat_filt_edge_coll(edge_t);

  Simplex_tree stree;
  mat_filt_edge_coll.filtered_edge_collapse(
    [&stree](std::vector<Vertex_handle> edge, Filtration_value filtration) {
        // insert the 2 vertices with a 0. filtration value just like a Rips
        stree.insert_simplex({edge[0]}, 0.);
        stree.insert_simplex({edge[1]}, 0.);
        // insert the edge
        stree.insert_simplex(edge, filtration);
      });

  stree.expansion(dim_max);

  std::cout << "The complex contains " << stree.num_simplices() << " simplices  after collapse. \n";
  std::cout << "   and has dimension " << stree.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  stree.initialize_filtration();
  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(stree);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(3);

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

void program_options(int argc, char* argv[], std::string& csv_matrix_file, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()(
      "input-file", po::value<std::string>(&csv_matrix_file),
      "Name of file containing a distance matrix. Can be square or lower triangular matrix. Separator is ';'.");

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
    std::cout << "of a Rips complex after edge collapse defined on a set of distance matrix.\n \n";
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
