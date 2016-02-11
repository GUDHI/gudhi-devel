/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA Sophia Antipolis-Méditerranée (France)
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

#include <gudhi/Simplex_tree.h>
#include <gudhi/Relaxed_witness_complex.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Persistent_cohomology.h>
#include "Landmark_choice_random_knn.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>
#include <algorithm>
#include <set>
#include <queue>
#include <iterator>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include "generators.h"
#include "output.h"

using namespace Gudhi;
using namespace Gudhi::witness_complex;
using namespace Gudhi::persistent_cohomology;

typedef std::vector<Point_d> Point_Vector;
typedef Relaxed_witness_complex< Simplex_tree<> > RelaxedWitnessComplex;

int main (int argc, char * const argv[])
{
  if (argc != 6) {
      std::cerr << "Usage: " << argv[0]
                << " nbP nbL dim alpha limD\n";
      return 0;
    }
  /*
  boost::filesystem::path p;
  for (; argc > 2; --argc, ++argv)
    p /= argv[1];
  */
  
  int nbP       = atoi(argv[1]);
  int nbL       = atoi(argv[2]);
  int dim       = atoi(argv[3]);
  double alpha  = atof(argv[4]);
  int limD      = atoi(argv[5]);
  //Construct the Simplex Tree
  clock_t start, end;

  // Construct the Simplex Tree
  Simplex_tree<> simplex_tree;

  // Read the point file
  Point_Vector point_vector;
  generate_points_sphere(point_vector, nbP, dim);
  std::cout << "Successfully generated " << point_vector.size() << " points.\n";
  std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

  // Choose landmarks
  std::vector<std::vector< int > > knn;
  std::vector<std::vector< FT > > distances;
  Gudhi::witness_complex::landmark_choice_by_random_knn(point_vector, nbL, alpha, limD, knn, distances);
  
  // Compute witness complex
  start = clock();
  RelaxedWitnessComplex(distances,
                        knn,
                        simplex_tree,
                        nbL,
                        alpha,
                        limD);
  end = clock();
  double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Witness complex for " << nbL << " landmarks took "
            << time << " s. \n";
  // std::cout << simplex_tree << "\n";
  
  // Compute the persistence diagram of the complex
  persistent_cohomology::Persistent_cohomology< Simplex_tree<>, Field_Zp > pcoh(simplex_tree, false);
  int p = 3;
  pcoh.init_coefficients( p ); //initilizes the coefficient field for homology
  pcoh.compute_persistent_cohomology( alpha/10 );
  pcoh.output_diagram();
}
