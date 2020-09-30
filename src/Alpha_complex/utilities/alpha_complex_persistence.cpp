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
#include <CGAL/Epeck_d.h>

#include <gudhi/Alpha_complex.h>
#include <gudhi/Persistent_cohomology.h>
// to construct a simplex_tree from alpha complex
#include <gudhi/Simplex_tree.h>
#include <gudhi/Points_off_io.h>

#include <iostream>
#include <string>
#include <limits>  // for numeric_limits
#include <vector>
#include <fstream>

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;

void program_options(int argc, char *argv[], std::string &off_file_points, bool &exact, bool &fast,
                     std::string &weight_file, std::string &output_file_diag, Filtration_value &alpha_square_max_value,
                     int &coeff_field_characteristic, Filtration_value &min_persistence);

template<class Point_d>
std::vector<Point_d> read_off(const std::string &off_file_points) {
  Gudhi::Points_off_reader<Point_d> off_reader(off_file_points);
  if (!off_reader.is_valid()) {
    std::cerr << "Alpha_complex - Unable to read file " << off_file_points << "\n";
    exit(-1);  // ----- >>
  }
  return off_reader.get_point_cloud();
}

std::vector<double> read_weight_file(const std::string &weight_file) {
  std::vector<double> weights;
  // Read weights information from file
  std::ifstream weights_ifstr(weight_file);
  if (weights_ifstr.good()) {
    double weight = 0.0;
    // Attempt read the weight in a double format, return false if it fails
    while (weights_ifstr >> weight) {
      weights.emplace_back(weight);
    }
  } else {
    std::cerr << "Unable to read weights file " << weight_file << std::endl;
    exit(-1);
  }
  return weights;
}

template<class Kernel>
Simplex_tree create_simplex_tree(const std::string &off_file_points, const std::string &weight_file,
                                 bool exact_version, Filtration_value alpha_square_max_value) {
  Simplex_tree stree;
  auto points = read_off<typename Kernel::Point_d>(off_file_points);

  if (weight_file != std::string()) {
    std::vector<double> weights = read_weight_file(weight_file);
    if (points.size() != weights.size()) {
      std::cerr << "Alpha_complex - Inconsistency between number of points (" << points.size()
                << ") and number of weights (" << weights.size() << ")" << "\n";
      exit(-1);  // ----- >>
    }
    // Init of an alpha complex from an OFF file
    Gudhi::alpha_complex::Alpha_complex<Kernel, true> alpha_complex_from_file(points, weights);

    if (!alpha_complex_from_file.create_complex(stree, alpha_square_max_value, exact_version)) {
      std::cerr << "Alpha complex simplicial complex creation failed." << std::endl;
      exit(-1);
    }
  } else {
    // Init of an alpha complex from an OFF file
    Gudhi::alpha_complex::Alpha_complex<Kernel> alpha_complex_from_file(points);

    if (!alpha_complex_from_file.create_complex(stree, alpha_square_max_value, exact_version)) {
      std::cerr << "Alpha complex simplicial complex creation failed." << std::endl;
      exit(-1);
    }
  }
  return stree;
}

int main(int argc, char **argv) {
  std::string weight_file;
  std::string off_file_points;
  std::string output_file_diag;
  bool exact_version = false;
  bool fast_version = false;
  Filtration_value alpha_square_max_value;
  int coeff_field_characteristic;
  Filtration_value min_persistence;

  program_options(argc, argv, off_file_points, exact_version, fast_version, weight_file, output_file_diag,
                  alpha_square_max_value, coeff_field_characteristic, min_persistence);

  if ((exact_version) && (fast_version)) {
    std::cerr << "You cannot set the exact and the fast version." << std::endl;
    exit(-1);
  }

  Simplex_tree stree;
  if (fast_version) {
    // WARNING : CGAL::Epick_d is fast but not safe (unlike CGAL::Epeck_d)
    // (i.e. when the points are on a grid)
    using Fast_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
    stree = create_simplex_tree<Fast_kernel>(off_file_points, weight_file, exact_version, alpha_square_max_value);
  } else {
    using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
    stree = create_simplex_tree<Kernel>(off_file_points, weight_file, exact_version, alpha_square_max_value);
  }
  // ----------------------------------------------------------------------------
  // Display information about the alpha complex
  // ----------------------------------------------------------------------------
  std::clog << "Simplicial complex is of dimension " << stree.dimension() << " - " << stree.num_simplices()
            << " simplices - " << stree.num_vertices() << " vertices." << std::endl;

  std::clog << "Simplex_tree dim: " << stree.dimension() << std::endl;
  // Compute the persistence diagram of the complex
  Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Gudhi::persistent_cohomology::Field_Zp> pcoh(
      stree);
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
                     std::string &weight_file, std::string &output_file_diag, Filtration_value &alpha_square_max_value,
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
      "weight-file,w", po::value<std::string>(&weight_file)->default_value(std::string()),
      "Name of file containing a point weights. Format is one weight per line:\n  W1\n  ...\n  Wn ")(
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
