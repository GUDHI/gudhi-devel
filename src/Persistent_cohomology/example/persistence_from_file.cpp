/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/reader_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>

#include <boost/program_options.hpp>

#include <string>

using namespace Gudhi;
using namespace Gudhi::persistent_cohomology;

typedef double Filtration_value;

void program_options(int argc, char * argv[]
                     , std::string & simplex_tree_file
                     , std::string & output_file
                     , int & p
                     , Filtration_value & min_persistence);

int main(int argc, char * argv[]) {
  std::string simplex_tree_file;
  std::string output_file;
  int p;
  Filtration_value min_persistence;

  program_options(argc, argv, simplex_tree_file, output_file, p, min_persistence);

  std::cout << "Simplex_tree from file=" << simplex_tree_file.c_str() << " - output_file=" << output_file.c_str()
      << std::endl;
  std::cout << "     - p=" << p << " - min_persistence=" << min_persistence << std::endl;

  // Read the list of simplices from a file.
  Simplex_tree<> simplex_tree;

  std::ifstream simplex_tree_stream(simplex_tree_file);
  simplex_tree_stream >> simplex_tree;

  std::cout << "The complex contains " << simplex_tree.num_simplices() << " simplices" << std::endl;
  std::cout << "   - dimension " << simplex_tree.dimension() << std::endl;

  /*
  std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for( auto f_simplex : simplex_tree.filtration_simplex_range() )
  { std::cout << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
  for( auto vertex : simplex_tree.simplex_vertex_range(f_simplex) )
  { std::cout << vertex << " "; }
  std::cout << std::endl;
  }*/

  // Sort the simplices in the order of the filtration
  simplex_tree.initialize_filtration();

  // Compute the persistence diagram of the complex
  Persistent_cohomology< Simplex_tree<>, Field_Zp > pcoh(simplex_tree);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(p);

  pcoh.compute_persistent_cohomology(min_persistence);

  // Output the diagram in output_file
  if (output_file.empty()) {
    pcoh.output_diagram();
  } else {
    std::ofstream out(output_file);
    pcoh.output_diagram(out);
    out.close();
  }

  return 0;
}

void program_options(int argc, char * argv[]
                     , std::string & simplex_tree_file
                     , std::string & output_file
                     , int & p
                     , Filtration_value & min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()
      ("input-file", po::value<std::string>(&simplex_tree_file),
       "Name of file containing a simplex set. Format is one simplex per line (cf. reader_utils.h - read_simplex): Dim1 X11 X12 ... X1d Fil1  ");

  po::options_description visible("Allowed options", 100);
  visible.add_options()
      ("help,h", "produce help message")
      ("output-file,o", po::value<std::string>(&output_file)->default_value(std::string()),
       "Name of file in which the persistence diagram is written. Default print in std::cout")
      ("field-charac,p", po::value<int>(&p)->default_value(11),
       "Characteristic p of the coefficient field Z/pZ for computing homology.")
      ("min-persistence,m", po::value<Filtration_value>(&min_persistence),
       "Minimal lifetime of homology feature to be recorded. Default is 0");

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
    exit(-1);
  }
}
