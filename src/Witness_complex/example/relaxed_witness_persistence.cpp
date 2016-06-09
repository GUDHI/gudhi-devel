/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016  INRIA Sophia Antipolis-Méditerranée (France)
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
#include <gudhi/Dim_lists.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Persistent_cohomology.h>
#include "Landmark_choice_random_knn.h"
#include "Landmark_choice_sparsification.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>
#include <algorithm>
#include <set>
#include <queue>
#include <iterator>
#include <string>

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
typedef Simplex_tree<>::Simplex_handle Simplex_handle;

/**
 * \brief Customized version of read_points
 * which takes into account a possible nbP first line
 *
 */
inline void
read_points_cust(std::string file_name, std::vector< std::vector< double > > & points) {
  std::ifstream in_file(file_name.c_str(), std::ios::in);
  if (!in_file.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return;
  }
  std::string line;
  double x;
  while (getline(in_file, line)) {
    std::vector< double > point;
    std::istringstream iss(line);
    while (iss >> x) {
      point.push_back(x);
    }
    if (point.size() != 1)
      points.push_back(point);
  }
  in_file.close();
}

void output_experiment_information(char * const file_name)
{
    std::cout << "Enter a valid experiment number. Usage: "
              << file_name << " exp_no options\n";
    std::cout << "Experiment description:\n"
              << "0 nbP nbL dim alpha limD mu_epsilon: "
              << "Build persistence diagram on relaxed witness complex "
              << "built from a point cloud on (dim-1)-dimensional sphere "
              << "consisting of nbP witnesses and nbL landmarks. "
              << "The maximal relaxation is alpha and the limit on simplicial complex "
              << "dimension is limD.\n";
    std::cout << "1 file_name nbL alpha limD: "
              << "Build persistence diagram on relaxed witness complex "
              << "build from a point cloud stored in a file and nbL landmarks. "
              << "The maximal relaxation is alpha and the limit on simplicial complex dimension is limD\n";
}

void rw_experiment(Point_Vector & point_vector, int nbL, FT alpha, int limD, FT mu_epsilon = 0.1)
{
  clock_t start, end;
  Simplex_tree<> simplex_tree;

  // Choose landmarks
  std::vector<std::vector< int > > knn;
  std::vector<std::vector< FT > > distances;
  start = clock();
  //Gudhi::witness_complex::landmark_choice_by_random_knn(point_vector, nbL, alpha, limD, knn, distances);
  std::vector<Point_d> landmarks;
  Gudhi::witness_complex::landmark_choice_by_sparsification(point_vector, nbL, mu_epsilon, landmarks);
  Gudhi::witness_complex::build_distance_matrix(point_vector,   // aka witnesses
                                                landmarks,  // aka landmarks
                                                alpha,
                                                limD,
                                                knn,
                                                distances);
  end = clock();
  double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Choice of " << nbL << " landmarks took "
            << time << " s. \n";
  // Compute witness complex
  start = clock();
  RelaxedWitnessComplex rw(distances,
                           knn,
                           simplex_tree,
                           nbL,
                           alpha,
                           limD);
  end = clock();
  time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Relaxed witness complex for " << nbL << " landmarks took "
            << time << " s. \n";
  std::cout << "The complex contains " << simplex_tree.num_simplices() << " simplices \n";
  // std::cout << simplex_tree << "\n";
  
  // Compute the persistence diagram of the complex
  simplex_tree.set_dimension(limD);
  persistent_cohomology::Persistent_cohomology< Simplex_tree<>, Field_Zp > pcoh(simplex_tree, true);
  int p = 3;
  pcoh.init_coefficients( p ); //initilizes the coefficient field for homology
  start = clock();
  pcoh.compute_persistent_cohomology( alpha/1000 );
  end = clock();
  time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Persistence diagram took "
            << time << " s. \n";
  pcoh.output_diagram();

  int chi = 0;
  for (auto sh: simplex_tree.complex_simplex_range())
    chi += 1-2*(simplex_tree.dimension(sh)%2);
  std::cout << "Euler characteristic is " << chi << std::endl;

}


int experiment0 (int argc, char * const argv[])
{
  if (argc != 8) {
    std::cerr << "Usage: " << argv[0]
              << " 0 nbP nbL dim alpha limD mu_epsilon\n";
    return 0;
  }
  /*
    boost::filesystem::path p;
    for (; argc > 2; --argc, ++argv)
    p /= argv[1];
  */
  
  int nbP       = atoi(argv[2]);
  int nbL       = atoi(argv[3]);
  int dim       = atoi(argv[4]);
  double alpha  = atof(argv[5]);
  int limD      = atoi(argv[6]);
  double mu_epsilon = atof(argv[7]);

  // Read the point file
  Point_Vector point_vector;
  generate_points_sphere(point_vector, nbP, dim);
  std::cout << "Successfully generated " << point_vector.size() << " points.\n";
  std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

  rw_experiment(point_vector, nbL, alpha, limD);
}

int experiment1 (int argc, char * const argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " 1 file_name\n";
    return 0;
  }
  /*
    boost::filesystem::path p;
    for (; argc > 2; --argc, ++argv)
    p /= argv[1];
  */
  
  std::string file_name = argv[2];

  // Read the point file
  Point_Vector point_vector;
  read_points_cust(file_name, point_vector);
  std::cout << "The file contains " << point_vector.size() << " points.\n";
  std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

  bool ok = false;
  int nbL, limD;
  double alpha;
  while (!ok) {
    std::cout << "Relaxed witness complex: parameters nbL, alpha, limD.\n";
    std::cout << "Enter nbL: ";
    std::cin >> nbL;
    std::cout << "Enter alpha: ";
    std::cin >> alpha;
    std::cout << "Enter limD: ";
    std::cin >> limD;
    std::cout << "Start relaxed witness complex...\n";
    rw_experiment(point_vector, nbL, alpha, limD);
    std::cout << "Is the result correct? [y/n]: ";
    char answer;
    std::cin >> answer;
    switch (answer) {
    case 'n':
      ok = false; break;
    default :
      ok = true; break;
    }
  }
  // ok = false;
  // while (!ok) {
  //   std::cout << "Rips complex: parameters threshold, limD.\n";
  //   std::cout << "Enter threshold: ";
  //   std::cin >> alpha;
  //   std::cout << "Enter limD: ";
  //   std::cin >> limD;
  //   std::cout << "Start Rips complex...\n";
  //   rips_experiment(point_vector, alpha, limD);
  //   std::cout << "Is the result correct? [y/n]: ";
  //   char answer;
  //   std::cin >> answer;
  //   switch (answer) {
  //   case 'n':
  //     ok = false; break;
  //   default :
  //     ok = true; break;
  //   }
  // }
}

int main (int argc, char * const argv[])
{
  if (argc == 1) {
    output_experiment_information(argv[0]);
    return 1;
  }
  switch (atoi(argv[1])) {
  case 0 :
    return experiment0(argc, argv);
  case 1 :
    return experiment1(argc, argv);
  default :
    output_experiment_information(argv[0]);
    return 1;
  }
}
