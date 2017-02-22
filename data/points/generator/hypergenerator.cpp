/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014  INRIA Saclay (France)
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

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>

#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>  // for std::ofstream
#include <cstdlib>

typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > K;
typedef K::Point_d Point;

void usage(char * const progName) {
  std::cerr << "Usage: " << progName << " in|on sphere|cube off_file_name points_number[integer > 0] " <<
      "dimension[integer > 1] radius[double > 0.0 | default = 1.0]" << std::endl;
  exit(-1);
}

int main(int argc, char **argv) {
  // program args management
  if ((argc != 6) && (argc != 7)) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct" << std::endl;
    usage(argv[0]);
  }

  int points_number = atoi(argv[4]);
  if (points_number <= 0) {
    std::cerr << "Error: " << argv[4] << " is not correct" << std::endl;
    usage(argv[0]);
  }

  int dimension = atoi(argv[5]);
  if (dimension <= 0) {
    std::cerr << "Error: " << argv[5] << " is not correct" << std::endl;
    usage(argv[0]);
  }

  double radius = 1.0;
  if (argc == 7) {
    radius = atof(argv[6]);
    if (radius <= 0.0) {
      std::cerr << "Error: " << argv[6] << " is not correct" << std::endl;
      usage(argv[0]);
    }
  }

  bool in = false;
  if (strcmp(argv[1], "in") == 0) {
    in = true;
  } else if (strcmp(argv[1], "on") != 0) {
      std::cerr << "Error: " << argv[1] << " is not correct" << std::endl;
      usage(argv[0]);
  }

  bool sphere = false;
  if (memcmp(argv[2], "sphere", sizeof("sphere")) == 0) {
    sphere = true;
  } else if (memcmp(argv[2], "cube", sizeof("cube")) != 0) {
      std::cerr << "Error: " << argv[2] << " is not correct" << std::endl;
      usage(argv[0]);
  }

  std::ofstream diagram_out(argv[3]);
  if (dimension == 3) {
    diagram_out << "OFF" << std::endl;
    diagram_out << points_number << " 0 0" << std::endl;
  } else {
    diagram_out << "nOFF" << std::endl;
    diagram_out << dimension << " " << points_number << " 0 0" << std::endl;
  }

  if (diagram_out.is_open()) {
    // Instanciate a random point generator
    CGAL::Random rng(0);
    // Generate "points_number" random points in a vector
    std::vector<Point> points;
    if (in) {
      if (sphere) {
        CGAL::Random_points_in_ball_d<Point> rand_it(dimension, radius, rng);
        CGAL::cpp11::copy_n(rand_it, points_number, std::back_inserter(points));
      } else {
        CGAL::Random_points_in_cube_d<Point> rand_it(dimension, radius, rng);
        CGAL::cpp11::copy_n(rand_it, points_number, std::back_inserter(points));
      }
    } else {  // means "on"
      if (sphere) {
        CGAL::Random_points_on_sphere_d<Point> rand_it(dimension, radius, rng);
        CGAL::cpp11::copy_n(rand_it, points_number, std::back_inserter(points));
      } else {
        std::cerr << "Sorry: on cube is not available" << std::endl;
        usage(argv[0]);
      }
    }

    for (auto thePoint : points) {
      int i = 0;
      for (; i < dimension - 1; i++) {
        diagram_out << thePoint[i] << " ";
      }
      diagram_out << thePoint[i] << std::endl;  // last point + Carriage Return
    }
  } else {
    std::cerr << "Error: " << argv[3] << " cannot be opened" << std::endl;
    usage(argv[0]);
  }

  return 0;
}

