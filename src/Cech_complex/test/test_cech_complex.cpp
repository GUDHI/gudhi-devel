/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018  Inria
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "cech_complex"
#include <boost/test/unit_test.hpp>

#include <cmath>  // float comparison
#include <limits>
#include <string>
#include <vector>
#include <algorithm>    // std::max

#include <gudhi/Cech_complex.h>
// to construct Cech_complex from a OFF file of points
#include <gudhi/Points_off_io.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Unitary_tests_utils.h>

#include <Miniball/Miniball.hpp>

// Type definitions
using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Point = std::vector<double>;
using Point_cloud = std::vector<Point>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;
using Cech_complex = Gudhi::cech_complex::Cech_complex<Simplex_tree, Point_cloud>;

using Point_iterator = Point_cloud::const_iterator;
using Coordinate_iterator = Point::const_iterator;
using Min_sphere = Miniball::Miniball<Miniball::CoordAccessor<Point_iterator, Coordinate_iterator>>;

BOOST_AUTO_TEST_CASE(Cech_complex_from_file) {
  // ----------------------------------------------------------------------------
  //
  // Init of a Cech complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("alphacomplexdoc.off");
  double threshold = 12.0;
  std::cout << "========== OFF FILE NAME = " << off_file_name << " - Cech threshold=" <<
      threshold << "==========" << std::endl;

  Points_off_reader off_reader(off_file_name);
  Point_cloud point_cloud = off_reader.get_point_cloud();
  Cech_complex cech_complex_from_file(point_cloud, threshold, Gudhi::Euclidean_distance());

  std::size_t i = 0;
  for (; i < point_cloud.size(); i++) {
    BOOST_CHECK(point_cloud[i] == *(cech_complex_from_file.point_iterator(i)));
  }
#ifdef GUDHI_DEBUG
  BOOST_CHECK_THROW (cech_complex_from_file.point_iterator(i+1), std::out_of_range);
#endif  // GUDHI_DEBUG

  const int DIMENSION_1 = 1;
  Simplex_tree st;
  cech_complex_from_file.create_complex(st, DIMENSION_1);
  std::cout << "st.dimension()=" << st.dimension() << std::endl;
  BOOST_CHECK(st.dimension() == DIMENSION_1);

  const int NUMBER_OF_VERTICES = 7;
  std::cout << "st.num_vertices()=" << st.num_vertices() << std::endl;
  BOOST_CHECK(st.num_vertices() == NUMBER_OF_VERTICES);

  std::cout << "st.num_simplices()=" << st.num_simplices() << std::endl;
  BOOST_CHECK(st.num_simplices() == 18);

  // Check filtration values of vertices is 0.0
  for (auto f_simplex : st.skeleton_simplex_range(0)) {
    BOOST_CHECK(st.filtration(f_simplex) == 0.0);
  }

  // Check filtration values of edges
  for (auto f_simplex : st.skeleton_simplex_range(DIMENSION_1)) {
    if (DIMENSION_1 == st.dimension(f_simplex)) {
      std::vector<Point> vp;
      std::cout << "vertex = (";
      for (auto vertex : st.simplex_vertex_range(f_simplex)) {
        std::cout << vertex << ",";
        vp.push_back(off_reader.get_point_cloud().at(vertex));
      }
      std::cout << ") - distance =" << Gudhi::Euclidean_distance()(vp.at(0), vp.at(1)) <<
          " - filtration =" << st.filtration(f_simplex) << std::endl;
      BOOST_CHECK(vp.size() == 2);
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex), Gudhi::Euclidean_distance()(vp.at(0), vp.at(1)));
    }
  }

  const int DIMENSION_2 = 2;
  Simplex_tree st2;
  cech_complex_from_file.create_complex(st2, DIMENSION_2);
  std::cout << "st2.dimension()=" << st2.dimension() << std::endl;
  BOOST_CHECK(st2.dimension() == DIMENSION_2);
  
  std::cout << "st2.num_vertices()=" << st2.num_vertices() << std::endl;
  BOOST_CHECK(st2.num_vertices() == NUMBER_OF_VERTICES);

  std::cout << "st2.num_simplices()=" << st2.num_simplices() << std::endl;
  BOOST_CHECK(st2.num_simplices() == 23);

  Point_cloud points012;
  for (std::size_t vertex = 0; vertex <= 2; vertex++) {
    points012.push_back(Point(cech_complex_from_file.point_iterator(vertex)->begin(),
                           cech_complex_from_file.point_iterator(vertex)->end()));
  }
  Min_sphere ms012(cech_complex_from_file.dimension(), points012.begin(),points012.end());

  Simplex_tree::Filtration_value f012 = st2.filtration(st2.find({0, 1, 2}));
  std::cout << "f012= " << f012 << " | ms012_radius= " << std::sqrt(ms012.squared_radius()) << std::endl;

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f012, std::sqrt(ms012.squared_radius()));

  Point_cloud points456;
  for (std::size_t vertex = 4; vertex <= 6; vertex++) {
    points456.push_back(Point(cech_complex_from_file.point_iterator(vertex)->begin(),
                           cech_complex_from_file.point_iterator(vertex)->end()));
  }
  Min_sphere ms456(cech_complex_from_file.dimension(), points456.begin(),points456.end());

  Simplex_tree::Filtration_value f456 = st2.filtration(st2.find({4, 5, 6}));
  std::cout << "f456= " << f456 << " | ms456_radius= " << std::sqrt(ms456.squared_radius()) << std::endl;

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f456, std::sqrt(ms456.squared_radius()));

  const int DIMENSION_3 = 3;
  Simplex_tree st3;
  cech_complex_from_file.create_complex(st3, DIMENSION_3);
  std::cout << "st3.dimension()=" << st3.dimension() << std::endl;
  BOOST_CHECK(st3.dimension() == DIMENSION_3);
  
  std::cout << "st3.num_vertices()=" << st3.num_vertices() << std::endl;
  BOOST_CHECK(st3.num_vertices() == NUMBER_OF_VERTICES);

  std::cout << "st3.num_simplices()=" << st3.num_simplices() << std::endl;
  BOOST_CHECK(st3.num_simplices() == 24);

  Point_cloud points0123;
  for (std::size_t vertex = 0; vertex <= 3; vertex++) {
    points0123.push_back(Point(cech_complex_from_file.point_iterator(vertex)->begin(),
                               cech_complex_from_file.point_iterator(vertex)->end()));
  }
  Min_sphere ms0123(cech_complex_from_file.dimension(), points0123.begin(),points0123.end());

  Simplex_tree::Filtration_value f0123 = st3.filtration(st3.find({0, 1, 2, 3}));
  std::cout << "f0123= " << f0123 << " | ms0123_radius= " << std::sqrt(ms0123.squared_radius()) << std::endl;

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f0123, std::sqrt(ms0123.squared_radius()));



  Point_cloud points01;
  for (std::size_t vertex = 0; vertex <= 1; vertex++) {
    points01.push_back(Point(cech_complex_from_file.point_iterator(vertex)->begin(),
                           cech_complex_from_file.point_iterator(vertex)->end()));
  }
  Min_sphere ms01(cech_complex_from_file.dimension(), points01.begin(),points01.end());

  Simplex_tree::Filtration_value f01 = st2.filtration(st2.find({0, 1}));
  std::cout << "f01= " << f01 << " | ms01_radius= " << std::sqrt(ms01.squared_radius()) << std::endl;

}

BOOST_AUTO_TEST_CASE(Cech_complex_from_points) {
  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Point_cloud points;
  std::vector<double> coords = { 0.0, 0.0, 0.0, 1.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 0.0, 0.0, 1.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 0.0, 1.0, 0.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 1.0, 0.0, 0.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));

  // ----------------------------------------------------------------------------
  // Init of a Cech complex from the list of points
  // ----------------------------------------------------------------------------
  Cech_complex cech_complex_from_points(points, 2.0, Gudhi::Euclidean_distance());

  std::cout << "========== cech_complex_from_points ==========" << std::endl;
  Simplex_tree st;
  const int DIMENSION = 3;
  cech_complex_from_points.create_complex(st, DIMENSION);

  // Another way to check num_simplices
  std::cout << "Iterator on Cech complex simplices in the filtration order, with [filtration value]:" << std::endl;
  int num_simplices = 0;
  for (auto f_simplex : st.filtration_simplex_range()) {
    num_simplices++;
    std::cout << "   ( ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << st.filtration(f_simplex) << "] ";
    std::cout << std::endl;
  }
  BOOST_CHECK(num_simplices == 15);
  std::cout << "st.num_simplices()=" << st.num_simplices() << std::endl;
  BOOST_CHECK(st.num_simplices() == 15);

  std::cout << "st.dimension()=" << st.dimension() << std::endl;
  BOOST_CHECK(st.dimension() == DIMENSION);
  std::cout << "st.num_vertices()=" << st.num_vertices() << std::endl;
  BOOST_CHECK(st.num_vertices() == 4);

  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "dimension(" << st.dimension(f_simplex) << ") - f = " << st.filtration(f_simplex) << std::endl;
    switch (st.dimension(f_simplex)) {
      case 0:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex), 0.0);
        break;
      case 1:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex), 1.41421, .00001);
        break;
      case 2:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex), 0.816497, .00001);
        break;
      case 3:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex), 0.866025, .00001);
        break;
      default:
        BOOST_CHECK(false);  // Shall not happen
        break;
    }
  }
}

#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE(Cech_create_complex_throw) {
  // ----------------------------------------------------------------------------
  //
  // Init of a Cech complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("alphacomplexdoc.off");
  double threshold = 12.0;
  std::cout << "========== OFF FILE NAME = " << off_file_name << " - Cech threshold=" <<
      threshold << "==========" << std::endl;

  Gudhi::Points_off_reader<Point> off_reader(off_file_name);
  Cech_complex cech_complex_from_file(off_reader.get_point_cloud(), threshold, Gudhi::Euclidean_distance());

  Simplex_tree stree;
  std::vector<int> simplex = {0, 1, 2};
  stree.insert_simplex_and_subfaces(simplex);
  std::cout << "Check exception throw in debug mode" << std::endl;
  // throw excpt because stree is not empty
  BOOST_CHECK_THROW (cech_complex_from_file.create_complex(stree, 1), std::invalid_argument);
}
#endif
