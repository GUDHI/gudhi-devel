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

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/Construct_closest_landmark_table.h>
#include <gudhi/Pick_random_points.h>
#include <gudhi/reader_utils.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>

typedef std::vector< Vertex_handle > typeVectorVertex;
typedef std::vector< std::vector <double> > Point_Vector;

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

int main(int argc, char * const argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_point_file nbL \n";
    return 0;
  }

  std::string file_name = argv[1];
  int nbL = atoi(argv[2]);
  clock_t start, end;

  // Construct the Simplex Tree
  Gudhi::Simplex_tree<> simplex_tree;

  // Read the point file
  Point_Vector point_vector, landmarks;
  read_points_cust(file_name, point_vector);
  std::cout << "Successfully read " << point_vector.size() << " points.\n";
  std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

  // Choose landmarks
  start = clock();
  std::vector<std::vector< int > > knn;
  Gudhi::pick_random_points(point_vector, 100, std::back_inserter(landmarks));
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
