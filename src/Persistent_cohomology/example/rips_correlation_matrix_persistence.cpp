/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko, Vincent Rouvreau
 *
 *    Copyright (C) 2016  INRIA
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

#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/reader_utils.h>

#include <boost/program_options.hpp>

#include <string>
#include <vector>
#include <limits>  // infinity


// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;
using Correlation_matrix = std::vector<std::vector<Filtration_value>>;

void program_options(int argc, char * argv[]
                     , std::string & csv_matrix_file
                     , std::string & filediag
                     , Filtration_value & threshold
                     , int & dim_max
                     , int & p
                     , Filtration_value & min_persistence);

int main(int argc, char * argv[]) {
  std::string csv_matrix_file;
  std::string filediag;
  Filtration_value threshold;
  int dim_max;
  int p;
  Filtration_value min_persistence;

  program_options(argc, argv, csv_matrix_file, filediag, threshold, dim_max, p, min_persistence);

  Correlation_matrix correlations = Gudhi::read_lower_triangular_matrix_from_csv_file<Filtration_value>(csv_matrix_file);
  
  //Given a correlation matrix M, we compute component-wise M'[i,j] = 1-M[i,j] to get a distance matrix:
  for ( size_t i = 0 ; i != correlations.size() ; ++i )
  {
	  for ( size_t j = 0 ; j != correlations[i].size() ; ++j )
	  {
		  correlations[i][j] = 1-correlations[i][j];
		  if ( correlations[i][j] < 0 )
		  {
			  std::cerr << "The input matrix is not a correlation matrix. \n";
			  throw "The input matrix is not a correlation matrix. \n";
		  }		  
	  }
  } 
  
  Rips_complex rips_complex_from_file(correlations, threshold);

  // Construct the Rips complex in a Simplex Tree
  Simplex_tree simplex_tree;

  rips_complex_from_file.create_complex(simplex_tree, dim_max);
  std::cout << "The complex contains " << simplex_tree.num_simplices() << " simplices \n";
  std::cout << "   and has dimension " << simplex_tree.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  simplex_tree.initialize_filtration();

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

void program_options(int argc, char * argv[]
                     , std::string & csv_matrix_file
                     , std::string & filediag
                     , Filtration_value & threshold
                     , int & dim_max
                     , int & p
                     , Filtration_value & min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()
      ("input-file", po::value<std::string>(&csv_matrix_file),
       "Name of file containing a distance matrix. Can be square or lower triangular matrix. Separator is ';'.");

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
    std::cout << "of a Rips complex defined on a set of distance matrix.\n \n";
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
