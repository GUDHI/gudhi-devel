/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/random_point_generators.h>

#include <CGAL/Epick_d.h>
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

  enum class Data_shape { sphere, cube, curve, torus, klein, undefined};

  Data_shape shape = Data_shape::undefined;
  if (memcmp(argv[2], "sphere", sizeof("sphere")) == 0) {
    shape = Data_shape::sphere;
  } else if (memcmp(argv[2], "cube", sizeof("cube")) == 0) {
    shape = Data_shape::cube;
  } else if (memcmp(argv[2], "curve", sizeof("curve")) == 0) {
    shape = Data_shape::curve;
  } else if (memcmp(argv[2], "torus", sizeof("torus")) == 0) {
    shape = Data_shape::torus;
  } else if (memcmp(argv[2], "klein", sizeof("klein")) == 0) {
    shape = Data_shape::klein;
  } else {
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
    // Generate "points_number" random points in a vector
    std::vector<Point> points;
    if (in) {
      switch (shape) {
        case Data_shape::sphere:
          points = Gudhi::generate_points_in_ball_d<K>(points_number, dimension, radius);
        break;
        case Data_shape::cube:
          points = Gudhi::generate_points_in_ball_d<K>(points_number, dimension, radius);
        break;
        case Data_shape::curve:
          std::cerr << "Sorry: in curve is not available" << std::endl;
          usage(argv[0]);
        break;
        case Data_shape::torus:
          std::cerr << "Sorry: in torus is not available" << std::endl;
          usage(argv[0]);
        break;
        case Data_shape::klein:
          std::cerr << "Sorry: in klein is not available" << std::endl;
          usage(argv[0]);
        break;
        default:
          usage(argv[0]);
        break;
      }
    } else {  // means "on"
      switch (shape) {
        case Data_shape::sphere:
          points = Gudhi::generate_points_on_sphere_d<K>(points_number, dimension, radius);
        break;
        case Data_shape::cube:
          std::cerr << "Sorry: on cube is not available" << std::endl;
          usage(argv[0]);
        break;
        case Data_shape::curve:
          points = Gudhi::generate_points_on_moment_curve<K>(points_number, dimension, -radius/2., radius/2.);
        break;
        case Data_shape::torus:
          if (dimension == 3)
            points = Gudhi::generate_points_on_torus_3D<K>(points_number, dimension, radius, radius/2.);
          else
            points = Gudhi::generate_points_on_torus_d<K>(points_number, dimension, true);
        break;
        case Data_shape::klein:
          switch (dimension) {
            case 3:
              points = Gudhi::generate_points_on_klein_bottle_3D<K>(points_number, radius, radius/2., true);
            break;
            case 4:
              points = Gudhi::generate_points_on_klein_bottle_4D<K>(points_number, radius, radius/2., 0., true);
            break;
            case 5:
              points = Gudhi::generate_points_on_klein_bottle_variant_5D<K>(points_number, radius, radius/2., true);
            break;
            default:
              std::cerr << "Sorry: on klein is only available for dimension 3, 4 and 5" << std::endl;
              usage(argv[0]);
            break;
          }
        break;
        default:
          usage(argv[0]);
        break;
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

