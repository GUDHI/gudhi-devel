/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2017 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Alpha_complex.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/Bottleneck.h>

#include <CGAL/Epick_d.h>

#include <boost/program_options.hpp>

#include <string>
#include <vector>
#include <limits>  // infinity
#include <utility>  // for pair
#include <algorithm>  // for transform


// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;
using Kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
using Point_d = Kernel::Point_d;
using Points_off_reader = Gudhi::Points_off_reader<Point_d>;

void program_options(int argc, char * argv[]
                     , std::string & off_file_points
                     , Filtration_value & threshold
                     , int & dim_max
                     , int & p
                     , Filtration_value & min_persistence);

static inline std::pair<double, double> compute_root_square(std::pair<double, double> input) {
  return std::make_pair(std::sqrt(input.first), std::sqrt(input.second));
}

int main(int argc, char * argv[]) {
  std::string off_file_points;
  Filtration_value threshold;
  int dim_max;
  int p;
  Filtration_value min_persistence;

  program_options(argc, argv, off_file_points, threshold, dim_max, p, min_persistence);

  Points_off_reader off_reader(off_file_points);

  // --------------------------------------------
  // Rips persistence
  // --------------------------------------------
  Rips_complex rips_complex(off_reader.get_point_cloud(), threshold, Gudhi::Euclidean_distance());

  // Construct the Rips complex in a Simplex Tree
  Simplex_tree rips_stree;

  rips_complex.create_complex(rips_stree, dim_max);
  std::clog << "The Rips complex contains " << rips_stree.num_simplices() << " simplices and has dimension "
            << rips_stree.dimension() << " \n";

  // Compute the persistence diagram of the complex
  Persistent_cohomology rips_pcoh(rips_stree);
  // initializes the coefficient field for homology
  rips_pcoh.init_coefficients(p);
  rips_pcoh.compute_persistent_cohomology(min_persistence);

  // rips_pcoh.output_diagram();

  // --------------------------------------------
  // Alpha persistence
  // --------------------------------------------
  Gudhi::alpha_complex::Alpha_complex<Kernel> alpha_complex(off_reader.get_point_cloud());

  Simplex_tree alpha_stree;
  alpha_complex.create_complex(alpha_stree, threshold * threshold);
  std::clog << "The Alpha complex contains " << alpha_stree.num_simplices() << " simplices and has dimension "
            << alpha_stree.dimension() << " \n";

  // Compute the persistence diagram of the complex
  Persistent_cohomology alpha_pcoh(alpha_stree);
  // initializes the coefficient field for homology
  alpha_pcoh.init_coefficients(p);
  alpha_pcoh.compute_persistent_cohomology(min_persistence * min_persistence);

  // alpha_pcoh.output_diagram();

  // --------------------------------------------
  // Bottleneck distance between both persistence
  // --------------------------------------------
  double max_b_distance {};
  for (int dim = 0; dim < dim_max; dim ++) {
    std::vector< std::pair< Filtration_value , Filtration_value > > rips_intervals;
    std::vector< std::pair< Filtration_value , Filtration_value > > alpha_intervals;
    rips_intervals = rips_pcoh.intervals_in_dimension(dim);
    alpha_intervals = alpha_pcoh.intervals_in_dimension(dim);
    std::transform(alpha_intervals.begin(), alpha_intervals.end(), alpha_intervals.begin(), compute_root_square);

    double bottleneck_distance = Gudhi::persistence_diagram::bottleneck_distance(rips_intervals, alpha_intervals);
    std::clog << "In dimension " << dim << ", bottleneck distance = " << bottleneck_distance << std::endl;
    if (bottleneck_distance > max_b_distance)
      max_b_distance = bottleneck_distance;
  }
  std::clog << "================================================================================" << std::endl;
  std::clog << "Bottleneck distance is " << max_b_distance << std::endl;

  return 0;
}

void program_options(int argc, char * argv[]
                     , std::string & off_file_points
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
