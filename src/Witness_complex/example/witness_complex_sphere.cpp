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
#define BOOST_PARAMETER_MAX_ARITY 12


#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/Construct_closest_landmark_table.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/reader_utils.h>

#include <CGAL/Epick_d.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>
#include <string>
#include <vector>

#include "generators.h"

/** Write a gnuplot readable file.
 *  Data range is a random access range of pairs (arg, value)
 */
template < typename Data_range >
void write_data(Data_range & data, std::string filename) {
  std::ofstream ofs(filename, std::ofstream::out);
  for (auto entry : data)
    ofs << entry.first << ", " << entry.second << "\n";
  ofs.close();
}

int main(int argc, char * const argv[]) {
  typedef Gudhi::witness_complex::Witness_complex<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>> Witness_complex;
  
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0]
        << " number_of_landmarks \n";
    return 0;
  }

  int number_of_landmarks = atoi(argv[1]);
  clock_t start, end;

  std::vector< std::pair<int, double> > l_time;

  // Read the point file
  for (int nbP = 500; nbP < 10000; nbP += 500) {
    // Construct the Simplex Tree
    Gudhi::Simplex_tree<> simplex_tree;
    Point_Vector point_vector, landmarks;
    generate_points_sphere(point_vector, nbP, 4);
    std::cout << "Successfully generated " << point_vector.size() << " points.\n";
    std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

    // Choose landmarks
    start = clock();
    Gudhi::subsampling::pick_n_random_points(point_vector, number_of_landmarks, std::back_inserter(landmarks));

    // Compute witness complex
    Witness_complex witness_complex(landmarks.begin(),
                                    landmarks.end(),
                                    point_vector.begin(),
                                    point_vector.end());
    witness_complex.create_complex(simplex_tree, 0);
    end = clock();
    double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Witness complex for " << number_of_landmarks << " landmarks took "
        << time << " s. \n";
    //assert(1 == 0);
    //std::cout << simplex_tree << "\n";
    std::cout << "Number of simplices is: " << simplex_tree.num_simplices() << "\n";
    l_time.push_back(std::make_pair(nbP, time));
  }
  write_data(l_time, "w_time.dat");
}
