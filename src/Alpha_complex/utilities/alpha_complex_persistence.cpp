/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <boost/program_options.hpp>

#include <CGAL/Epick_d.h>

#include <gudhi/Alpha_complex.h>
#include <gudhi/Persistent_cohomology.h>
// to construct a simplex_tree from alpha complex
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <string>
#include <limits>  // for numeric_limits

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;

void program_options(int argc, char *argv[], std::string &off_file_points, bool &exact, bool &fast,
                     std::string &output_file_diag, Filtration_value &alpha_square_max_value,
                     int &coeff_field_characteristic, Filtration_value &min_persistence);

int main(int argc, char **argv) {
  std::string off_file_points;
  std::string output_file_diag;
  bool exact_version = false;
  bool fast_version = false;
  Filtration_value alpha_square_max_value;
  int coeff_field_characteristic;
  Filtration_value min_persistence;

  program_options(argc, argv, off_file_points, exact_version, fast_version, output_file_diag, alpha_square_max_value,
                  coeff_field_characteristic, min_persistence);

  if ((exact_version) && (fast_version)) {
    std::cerr << "You cannot set the exact and the fast version." << std::endl;
    exit(-1);
  }

  Simplex_tree simplex;
  if (fast_version) {
    // WARNING : CGAL::Epick_d is fast but not safe (unlike CGAL::Epeck_d)
    // (i.e. when the points are on a grid)
    using Fast_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;

    // Init of an alpha complex from an OFF file
    Gudhi::alpha_complex::Alpha_complex<Fast_kernel> alpha_complex_from_file(off_file_points);

    if (!alpha_complex_from_file.create_complex(simplex, alpha_square_max_value)) {
      std::cerr << "Fast Alpha complex simplicial complex creation failed." << std::endl;
      exit(-1);
    }
  } else {
    using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;

    // Init of an alpha complex from an OFF file
    Gudhi::alpha_complex::Alpha_complex<Kernel> alpha_complex_from_file(off_file_points);

    if (!alpha_complex_from_file.create_complex(simplex, alpha_square_max_value, exact_version)) {
      std::cerr << "Alpha complex simplicial complex creation failed." << std::endl;
      exit(-1);
    }
  }
  // ----------------------------------------------------------------------------
  // Display information about the alpha complex
  // ----------------------------------------------------------------------------
  std::clog << "Simplicial complex is of dimension " << simplex.dimension() << " - " << simplex.num_simplices()
            << " simplices - " << simplex.num_vertices() << " vertices." << std::endl;

  std::clog << "Simplex_tree dim: " << simplex.dimension() << std::endl;
  // Compute the persistence diagram of the complex
  Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Gudhi::persistent_cohomology::Field_Zp> pcoh(
      simplex);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(coeff_field_characteristic);

  pcoh.compute_persistent_cohomology(min_persistence);

  // Output the diagram in filediag
  if (output_file_diag.empty()) {
    pcoh.output_diagram();
  } else {
    std::clog << "Result in file: " << output_file_diag << std::endl;
    std::ofstream out(output_file_diag);
    pcoh.output_diagram(out);
    out.close();
  }
  return 0;
}

void program_options(int argc, char *argv[], std::string &off_file_points, bool &exact, bool &fast,
                     std::string &output_file_diag, Filtration_value &alpha_square_max_value,
                     int &coeff_field_characteristic, Filtration_value &min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
                       "Name of file containing a point set. Format is one point per line:   X1 ... Xd ");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")(
      "exact,e", po::bool_switch(&exact),
      "To activate exact version of Alpha complex (default is false, not available if fast is set)")(
      "fast,f", po::bool_switch(&fast),
      "To activate fast version of Alpha complex (default is false, not available if exact is set)")(
      "output-file,o", po::value<std::string>(&output_file_diag)->default_value(std::string()),
      "Name of file in which the persistence diagram is written. Default print in std::clog")(
      "max-alpha-square-value,r", po::value<Filtration_value>(&alpha_square_max_value)
                                      ->default_value(std::numeric_limits<Filtration_value>::infinity()),
      "Maximal alpha square value for the Alpha complex construction.")(
      "field-charac,p", po::value<int>(&coeff_field_characteristic)->default_value(11),
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
    std::clog << "of an Alpha complex defined on a set of input points.\n \n";
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
