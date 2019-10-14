/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
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

#include <boost/program_options.hpp>
#include <boost/variant.hpp>

#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_3D_off_io.h>

#include <fstream>
#include <string>
#include <vector>
#include <limits>  // for numeric_limits<>

// gudhi type definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = Simplex_tree::Filtration_value;
using Persistent_cohomology =
    Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Gudhi::persistent_cohomology::Field_Zp>;

void program_options(int argc, char *argv[], std::string &off_file_points, bool &exact, bool &safe,
                     std::string &weight_file, std::string &cuboid_file, std::string &output_file_diag,
                     Filtration_value &alpha_square_max_value, int &coeff_field_characteristic,
                     Filtration_value &min_persistence);

bool read_weight_file(const std::string &weight_file, std::vector<double> &weights) {
  // Read weights information from file
  std::ifstream weights_ifstr(weight_file);
  if (weights_ifstr.good()) {
    double weight = 0.0;
    // Attempt read the weight in a double format, return false if it fails
    while (weights_ifstr >> weight) {
      weights.push_back(weight);
    }
  } else {
    return false;
  }
  return true;
}

bool read_cuboid_file(const std::string &cuboid_file, double &x_min, double &y_min, double &z_min, double &x_max,
                      double &y_max, double &z_max) {
  // Read weights information from file
  std::ifstream iso_cuboid_str(cuboid_file);
  if (iso_cuboid_str.is_open()) {
    if (!(iso_cuboid_str >> x_min >> y_min >> z_min >> x_max >> y_max >> z_max)) {
      return false;
    }
  } else {
    return false;
  }
  return true;
}

template <typename AlphaComplex3d>
std::vector<typename AlphaComplex3d::Point_3> read_off(const std::string &off_file_points) {
  // Read the OFF file (input file name given as parameter) and triangulate points
  Gudhi::Points_3D_off_reader<typename AlphaComplex3d::Point_3> off_reader(off_file_points);
  // Check the read operation was correct
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read OFF file " << off_file_points << std::endl;
    exit(-1);
  }
  return off_reader.get_point_cloud();
}

int main(int argc, char **argv) {
  std::string off_file_points;
  std::string weight_file;
  std::string cuboid_file;
  std::string output_file_diag;
  Filtration_value alpha_square_max_value = 0.;
  int coeff_field_characteristic = 0;
  Filtration_value min_persistence = 0.;
  bool exact_version = false;
  bool fast_version = false;
  bool weighted_version = false;
  bool periodic_version = false;

  program_options(argc, argv, off_file_points, exact_version, fast_version, weight_file, cuboid_file, output_file_diag,
                  alpha_square_max_value, coeff_field_characteristic, min_persistence);

  std::vector<double> weights;
  if (weight_file != std::string()) {
    if (!read_weight_file(weight_file, weights)) {
      std::cerr << "Unable to read weights file " << weight_file << std::endl;
      exit(-1);
    }
    weighted_version = true;
  }

  double x_min = 0., y_min = 0., z_min = 0., x_max = 0., y_max = 0., z_max = 0.;
  std::ifstream iso_cuboid_str(argv[3]);
  if (cuboid_file != std::string()) {
    if (!read_cuboid_file(cuboid_file, x_min, y_min, z_min, x_max, y_max, z_max)) {
      std::cerr << "Unable to read cuboid file " << cuboid_file << std::endl;
      exit(-1);
    }
    periodic_version = true;
  }

  Gudhi::alpha_complex::complexity complexity = Gudhi::alpha_complex::complexity::SAFE;
  if (exact_version) {
    if (fast_version) {
      std::cerr << "You cannot set the exact and the fast version." << std::endl;
      exit(-1);
    }
    complexity = Gudhi::alpha_complex::complexity::EXACT;
  }
  if (fast_version) {
    complexity = Gudhi::alpha_complex::complexity::FAST;
  }

  Simplex_tree simplex_tree;

  switch (complexity) {
    case Gudhi::alpha_complex::complexity::FAST:
      if (weighted_version) {
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, true>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points, weights, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, false>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points, weights);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      } else {
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, true>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, false>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      }
      break;
    case Gudhi::alpha_complex::complexity::EXACT:
      if (weighted_version) {
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, true>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points, weights, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, false>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points, weights);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      } else {
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, true>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, false>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      }
      break;
    case Gudhi::alpha_complex::complexity::SAFE:
      if (weighted_version) {
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, true>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points, weights, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, false>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points, weights);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      } else {
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, true>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, false>;
          auto points = read_off<Alpha_complex_3d>(off_file_points);
          Alpha_complex_3d alpha_complex(points);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      }
      break;
    default:
      std::cerr << "Unknown complexity value " << std::endl;
      exit(-1);
      break;
  }

  // Sort the simplices in the order of the filtration
  simplex_tree.initialize_filtration();

  std::cout << "Simplex_tree dim: " << simplex_tree.dimension() << std::endl;
  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(simplex_tree, true);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(coeff_field_characteristic);

  pcoh.compute_persistent_cohomology(min_persistence);

  // Output the diagram in filediag
  if (output_file_diag.empty()) {
    pcoh.output_diagram();
  } else {
    std::cout << "Result in file: " << output_file_diag << std::endl;
    std::ofstream out(output_file_diag);
    pcoh.output_diagram(out);
    out.close();
  }

  return 0;
}

void program_options(int argc, char *argv[], std::string &off_file_points, bool &exact, bool &fast,
                     std::string &weight_file, std::string &cuboid_file, std::string &output_file_diag,
                     Filtration_value &alpha_square_max_value, int &coeff_field_characteristic,
                     Filtration_value &min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
                       "Name of file containing a point set. Format is one point per line:   X1 ... Xd ");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")(
      "exact,e", po::bool_switch(&exact),
      "To activate exact version of Alpha complex 3d (default is false, not available if fast is set)")(
      "fast,f", po::bool_switch(&fast),
      "To activate fast version of Alpha complex 3d (default is false, not available if exact is set)")(
      "weight-file,w", po::value<std::string>(&weight_file)->default_value(std::string()),
      "Name of file containing a point weights. Format is one weight per line:\n  W1\n  ...\n  Wn ")(
      "cuboid-file,c", po::value<std::string>(&cuboid_file),
      "Name of file describing the periodic domain. Format is:\n  min_hx min_hy min_hz\n  max_hx max_hy max_hz")(
      "output-file,o", po::value<std::string>(&output_file_diag)->default_value(std::string()),
      "Name of file in which the persistence diagram is written. Default print in std::cout")(
      "max-alpha-square-value,r",
      po::value<Filtration_value>(&alpha_square_max_value)
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

  if (vm.count("help") || !vm.count("input-file") || !vm.count("weight-file")) {
    std::cout << std::endl;
    std::cout << "Compute the persistent homology with coefficient field Z/pZ \n";
    std::cout << "of a 3D Alpha complex defined on a set of input points.\n";
    std::cout << "3D Alpha complex can be safe (by default) exact or fast, weighted and/or periodic\n\n";
    std::cout << "The output diagram contains one bar per line, written with the convention: \n";
    std::cout << "   p   dim b d \n";
    std::cout << "where dim is the dimension of the homological feature,\n";
    std::cout << "b and d are respectively the birth and death of the feature and \n";
    std::cout << "p is the characteristic of the field Z/pZ used for homology coefficients.\n\n";

    std::cout << "Usage: " << argv[0] << " [options] input-file weight-file\n\n";
    std::cout << visible << std::endl;
    exit(-1);
  }
}
