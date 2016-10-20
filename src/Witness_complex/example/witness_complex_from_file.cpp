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

#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Points_off_io.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/Construct_closest_landmark_table.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/reader_utils.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>

typedef std::vector< Vertex_handle > typeVectorVertex;
typedef std::vector< std::vector <double> > Point_Vector;

int main(int argc, char * const argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_point_file.off nbL \n";
    return 0;
  }

  std::string off_file_name = argv[1];
  int nbL = atoi(argv[2]);
  clock_t start, end;

  // Construct the Simplex Tree
  Gudhi::Simplex_tree<> simplex_tree;

  // Read the OFF file (input file name given as parameter) and triangulate points
  Gudhi::Points_off_reader<std::vector <double>> off_reader(off_file_name);
  // Check the read operation was correct
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << off_file_name << std::endl;
  }
  // Read the point file
  Point_Vector point_vector = off_reader.get_point_cloud();
  std::cout << "Successfully read " << point_vector.size() << " points.\n";
  std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

  // Choose landmarks
  start = clock();
  std::vector<std::vector< int > > knn;
  Point_Vector landmarks;
  Gudhi::subsampling::pick_n_random_points(point_vector, 100, std::back_inserter(landmarks));
  Gudhi::witness_complex::construct_closest_landmark_table(point_vector, landmarks, knn);
  end = clock();
  std::cout << "Landmark choice for " << nbL << " landmarks took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";

  // Compute witness complex
  start = clock();
  Gudhi::witness_complex::witness_complex(knn, nbL, point_vector[0].size(), simplex_tree);
  end = clock();
  std::cout << "Witness complex took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
}
