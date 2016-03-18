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

#ifndef EXAMPLE_WITNESS_COMPLEX_GENERATORS_H_
#define EXAMPLE_WITNESS_COMPLEX_GENERATORS_H_

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>

#include <fstream>
#include <string>
#include <vector>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::FT FT;
typedef K::Point_d Point_d;
typedef std::vector<Point_d> Point_Vector;
typedef CGAL::Random_points_in_cube_d<Point_d> Random_cube_iterator;
typedef CGAL::Random_points_in_ball_d<Point_d> Random_point_iterator;

/**
 * \brief Rock age method of reading off file
 *
 */
inline void
off_reader_cust(std::string file_name, std::vector<Point_d> & points) {
  std::ifstream in_file(file_name.c_str(), std::ios::in);
  if (!in_file.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return;
  }
  std::string line;
  double x;
  // Line OFF. No need in it
  if (!getline(in_file, line)) {
    std::cerr << "No line OFF\n";
    return;
  }
  // Line with 3 numbers. No need
  if (!getline(in_file, line)) {
    std::cerr << "No line with 3 numbers\n";
    return;
  }
  // Reading points
  while (getline(in_file, line)) {
    std::vector< double > point;
    std::istringstream iss(line);
    while (iss >> x) {
      point.push_back(x);
    }
    points.push_back(Point_d(point));
  }
  in_file.close();
}

/**
 * \brief Customized version of read_points
 * which takes into account a possible nbP first line
 *
 */
inline void
read_points_cust(std::string file_name, Point_Vector & points) {
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
    Point_d p(point.begin(), point.end());
    if (point.size() != 1)
      points.push_back(p);
  }
  in_file.close();
}

/** \brief Generate points on a grid in a cube of side 2
 * having {+-1}^D as vertices and insert them in W.
 * The grid has "width" points on each side. 
 * If torus is true then it is supposed that the cube represents
 * a flat torus, hence the opposite borders are associated.
 * The points on border in this case are not placed twice.
 */
void generate_points_grid(Point_Vector& W, int width, int D, bool torus) {
  int nb_points = 1;
  for (int i = 0; i < D; ++i)
    nb_points *= width;
  for (int i = 0; i < nb_points; ++i) {
    std::vector<double> point;
    int cell_i = i;
    for (int l = 0; l < D; ++l) {
      if (torus)
        point.push_back(-1 + (2.0 / (width - 1))*(cell_i % width));
      else
        point.push_back(-1 + (2.0 / width)*(cell_i % width));
      // attention: the bottom and the right are covered too!
      cell_i /= width;
    }
    W.push_back(point);
  }
}

/** \brief Generate nbP points uniformly in a cube of side 2
 * having {+-1}^dim as its vertices and insert them in W.
 */
void generate_points_random_box(Point_Vector& W, int nbP, int dim) {
  Random_cube_iterator rp(dim, 1.0);
  for (int i = 0; i < nbP; i++) {
    W.push_back(*rp++);
  }
}

/** \brief Generate nbP points uniformly on a (dim-1)-sphere
 *  and insert them in W.
 */
void generate_points_sphere(Point_Vector& W, int nbP, int dim) {
  CGAL::Random_points_on_sphere_d<Point_d> rp(dim, 1);
  for (int i = 0; i < nbP; i++)
    W.push_back(*rp++);
}

#endif  // EXAMPLE_WITNESS_COMPLEX_GENERATORS_H_
