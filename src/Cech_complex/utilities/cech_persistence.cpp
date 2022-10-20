/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Cech_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_off_io.h>

#include <boost/program_options.hpp>

#include <CGAL/Epeck_d.h>  // For EXACT or SAFE version
#include <CGAL/Epick_d.h>  // For FAST version

#include <string>
#include <vector>
#include <limits>  // infinity

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = Simplex_tree::Filtration_value;

using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

void program_options(int argc, char* argv[], std::string& off_file_points, bool& exact, bool& fast,
                     std::string& filediag, Filtration_value& max_radius, int& dim_max, int& p,
                     Filtration_value& min_persistence);

template<class Kernel>
Simplex_tree create_simplex_tree(const std::string &off_file_points, bool exact_version,
                                 Filtration_value max_radius, int dim_max) {
  using Point = typename Kernel::Point_d;
  using Points_off_reader = Gudhi::Points_off_reader<Point>;
  using Cech_complex = Gudhi::cech_complex::Cech_complex<Kernel, Simplex_tree>;

  Simplex_tree stree;

  Points_off_reader off_reader(off_file_points);
  Cech_complex cech_complex_from_file(off_reader.get_point_cloud(), max_radius, exact_version);
  cech_complex_from_file.create_complex(stree, dim_max);

  return stree;
}

int main(int argc, char* argv[]) {
  std::string off_file_points;
  std::string filediag;
  bool exact_version = false;
  bool fast_version = false;
  Filtration_value max_radius;
  int dim_max;
  int p;
  Filtration_value min_persistence;

  program_options(argc, argv, off_file_points, exact_version, fast_version, filediag, max_radius, dim_max, p,
                  min_persistence);

  if ((exact_version) && (fast_version)) {
    std::cerr << "You cannot set the exact and the fast version." << std::endl;
    exit(-1);
  }

  Simplex_tree stree;
  if (fast_version) {
    // WARNING : CGAL::Epick_d is fast but not safe (unlike CGAL::Epeck_d)
    // (i.e. when the points are on a grid)
    using Fast_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
    stree = create_simplex_tree<Fast_kernel>(off_file_points, exact_version, max_radius, dim_max);
  } else {
    using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
    stree = create_simplex_tree<Kernel>(off_file_points, exact_version, max_radius, dim_max);
  }

  std::clog << "The complex contains " << stree.num_simplices() << " simplices \n";
  std::clog << "   and has dimension " << stree.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  stree.initialize_filtration();

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(stree);
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

void program_options(int argc, char* argv[], std::string& off_file_points, bool& exact, bool& fast,
                     std::string& filediag, Filtration_value& max_radius, int& dim_max, int& p,
                     Filtration_value& min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
                       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")(
      "exact,e", po::bool_switch(&exact),
      "To activate exact version of Cech complex (default is false, not available if fast is set)")(
      "fast,f", po::bool_switch(&fast),
      "To activate fast version of Cech complex (default is false, not available if exact is set)")(
      "output-file,o", po::value<std::string>(&filediag)->default_value(std::string()),
      "Name of file in which the persistence diagram is written. Default print in standard output")(
      "max-radius,r",
      po::value<Filtration_value>(&max_radius)->default_value(std::numeric_limits<Filtration_value>::infinity()),
      "Maximal length of an edge for the Cech complex construction.")(
      "cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
      "Maximal dimension of the Cech complex we want to compute.")(
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
    std::clog << "of a Cech complex defined on a set of input points.\n \n";
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
