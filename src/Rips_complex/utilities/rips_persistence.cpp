/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_off_io.h>

#include <boost/program_options.hpp>

#include <string>
#include <vector>
#include <limits>  // infinity

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
using Point = std::vector<double>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence);

int main(int argc, char* argv[]) {
  std::string off_file_points;
  std::string filediag;
  Filtration_value threshold;
  int dim_max;
  int p;
  Filtration_value min_persistence;

  program_options(argc, argv, off_file_points, filediag, threshold, dim_max, p, min_persistence);

  Points_off_reader off_reader(off_file_points);
  Rips_complex rips_complex_from_file(off_reader.get_point_cloud(), threshold, Gudhi::Euclidean_distance());

  // Construct the Rips complex in a Simplex Tree
  Simplex_tree simplex_tree;

  rips_complex_from_file.create_complex(simplex_tree, dim_max);
  std::clog << "The complex contains " << simplex_tree.num_simplices() << " simplices \n";
  std::clog << "   and has dimension " << simplex_tree.dimension() << " \n";

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(simplex_tree);
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

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
                       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")(
      "output-file,o", po::value<std::string>(&filediag)->default_value(std::string()),
      "Name of file in which the persistence diagram is written. Default print in std::clog")(
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
    std::clog << std::endl;
    std::clog << "Compute the persistent homology with coefficient field Z/pZ \n";
    std::clog << "of a Rips complex defined on a set of input points.\n \n";
    std::clog << "The output diagram contains one bar per line, written with the convention: \n";
    std::clog << "   p   dim b d \n";
    std::clog << "where dim is the dimension of the homological feature,\n";
    std::clog << "b and d are respectively the birth and death of the feature and \n";
    std::clog << "p is the characteristic of the field Z/pZ used for homology coefficients." << std::endl << std::endl;

    std::clog << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::clog << visible << std::endl;
    exit(-1);
  }
}
