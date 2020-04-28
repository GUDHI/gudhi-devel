/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "cech_complex"
#include <boost/test/unit_test.hpp>

#include <cmath>  // float comparison
#include <limits>
#include <string>
#include <vector>
#include <algorithm>  // std::max

#include <gudhi/Cech_complex.h>
// to construct Cech_complex from a OFF file of points
#include <gudhi/Points_off_io.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Unitary_tests_utils.h>
#include <gudhi/Miniball.hpp>

// Type definitions
using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Point = std::vector<Filtration_value>;
using Point_cloud = std::vector<Point>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;
using Cech_complex = Gudhi::cech_complex::Cech_complex<Simplex_tree, Point_cloud>;

using Point_iterator = Point_cloud::const_iterator;
using Coordinate_iterator = Point::const_iterator;
using Min_sphere = Gudhi::Miniball::Miniball<Gudhi::Miniball::CoordAccessor<Point_iterator, Coordinate_iterator>>;

BOOST_AUTO_TEST_CASE(Cech_complex_for_documentation) {
  // ----------------------------------------------------------------------------
  //
  // Init of a Cech complex from a point cloud
  //
  // ----------------------------------------------------------------------------
  Point_cloud points;
  points.push_back({1., 0.});                  // 0
  points.push_back({0., 1.});                  // 1
  points.push_back({2., 1.});                  // 2
  points.push_back({3., 2.});                  // 3
  points.push_back({0., 3.});                  // 4
  points.push_back({3. + std::sqrt(3.), 3.});  // 5
  points.push_back({1., 4.});                  // 6
  points.push_back({3., 4.});                  // 7
  points.push_back({2., 4. + std::sqrt(3.)});  // 8
  points.push_back({0., 4.});                  // 9
  points.push_back({-0.5, 2.});                // 10

  Filtration_value max_radius = 1.0;
  std::clog << "========== NUMBER OF POINTS = " << points.size() << " - Cech max_radius = " << max_radius
            << "==========" << std::endl;

  Cech_complex cech_complex_for_doc(points, max_radius);

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(cech_complex_for_doc.max_radius(), max_radius);
  std::size_t i = 0;
  for (; i < points.size(); i++) {
    BOOST_CHECK(points[i] == cech_complex_for_doc.get_point(i));
  }

  const int DIMENSION_1 = 1;
  Simplex_tree st;
  cech_complex_for_doc.create_complex(st, DIMENSION_1);
  std::clog << "st.dimension()=" << st.dimension() << std::endl;
  BOOST_CHECK(st.dimension() == DIMENSION_1);

  const int NUMBER_OF_VERTICES = 11;
  std::clog << "st.num_vertices()=" << st.num_vertices() << std::endl;
  BOOST_CHECK(st.num_vertices() == NUMBER_OF_VERTICES);

  std::clog << "st.num_simplices()=" << st.num_simplices() << std::endl;
  BOOST_CHECK(st.num_simplices() == 27);

  // Check filtration values of vertices is 0.0
  for (auto f_simplex : st.skeleton_simplex_range(0)) {
    BOOST_CHECK(st.filtration(f_simplex) == 0.0);
  }

  // Check filtration values of edges
  for (auto f_simplex : st.skeleton_simplex_range(DIMENSION_1)) {
    if (DIMENSION_1 == st.dimension(f_simplex)) {
      std::vector<Point> vp;
      std::clog << "vertex = (";
      for (auto vertex : st.simplex_vertex_range(f_simplex)) {
        std::clog << vertex << ",";
        vp.push_back(points.at(vertex));
      }
      std::clog << ") - distance =" << Gudhi::Minimal_enclosing_ball_radius()(vp.at(0), vp.at(1))
                << " - filtration =" << st.filtration(f_simplex) << std::endl;
      BOOST_CHECK(vp.size() == 2);
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex),
                                      Gudhi::Minimal_enclosing_ball_radius()(vp.at(0), vp.at(1)));
    }
  }

  const int DIMENSION_2 = 2;

#ifdef GUDHI_DEBUG
  BOOST_CHECK_THROW(cech_complex_for_doc.create_complex(st, DIMENSION_2), std::invalid_argument);
#endif

  Simplex_tree st2;
  cech_complex_for_doc.create_complex(st2, DIMENSION_2);
  std::clog << "st2.dimension()=" << st2.dimension() << std::endl;
  BOOST_CHECK(st2.dimension() == DIMENSION_2);

  std::clog << "st2.num_vertices()=" << st2.num_vertices() << std::endl;
  BOOST_CHECK(st2.num_vertices() == NUMBER_OF_VERTICES);

  std::clog << "st2.num_simplices()=" << st2.num_simplices() << std::endl;
  BOOST_CHECK(st2.num_simplices() == 30);

  Point_cloud points012;
  for (std::size_t vertex = 0; vertex <= 2; vertex++) {
    points012.push_back(cech_complex_for_doc.get_point(vertex));
  }
  std::size_t dimension = points[0].end() - points[0].begin();
  Min_sphere ms012(dimension, points012.begin(), points012.end());

  Simplex_tree::Filtration_value f012 = st2.filtration(st2.find({0, 1, 2}));
  std::clog << "f012= " << f012 << " | ms012_radius= " << std::sqrt(ms012.squared_radius()) << std::endl;

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f012, std::sqrt(ms012.squared_radius()));

  Point_cloud points1410;
  points1410.push_back(cech_complex_for_doc.get_point(1));
  points1410.push_back(cech_complex_for_doc.get_point(4));
  points1410.push_back(cech_complex_for_doc.get_point(10));
  Min_sphere ms1410(dimension, points1410.begin(), points1410.end());

  Simplex_tree::Filtration_value f1410 = st2.filtration(st2.find({1, 4, 10}));
  std::clog << "f1410= " << f1410 << " | ms1410_radius= " << std::sqrt(ms1410.squared_radius()) << std::endl;

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f1410, std::sqrt(ms1410.squared_radius()));

  Point_cloud points469;
  points469.push_back(cech_complex_for_doc.get_point(4));
  points469.push_back(cech_complex_for_doc.get_point(6));
  points469.push_back(cech_complex_for_doc.get_point(9));
  Min_sphere ms469(dimension, points469.begin(), points469.end());

  Simplex_tree::Filtration_value f469 = st2.filtration(st2.find({4, 6, 9}));
  std::clog << "f469= " << f469 << " | ms469_radius= " << std::sqrt(ms469.squared_radius()) << std::endl;

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f469, std::sqrt(ms469.squared_radius()));

  BOOST_CHECK((st2.find({6, 7, 8}) == st2.null_simplex()));
  BOOST_CHECK((st2.find({3, 5, 7}) == st2.null_simplex()));
}

BOOST_AUTO_TEST_CASE(Cech_complex_from_points) {
  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Point_cloud points;
  std::vector<double> coords = {0.0, 0.0, 0.0, 1.0};
  points.push_back(Point(coords.begin(), coords.end()));
  coords = {0.0, 0.0, 1.0, 0.0};
  points.push_back(Point(coords.begin(), coords.end()));
  coords = {0.0, 1.0, 0.0, 0.0};
  points.push_back(Point(coords.begin(), coords.end()));
  coords = {1.0, 0.0, 0.0, 0.0};
  points.push_back(Point(coords.begin(), coords.end()));

  // ----------------------------------------------------------------------------
  // Init of a Cech complex from the list of points
  // ----------------------------------------------------------------------------
  Cech_complex cech_complex_from_points(points, 2.0);

  std::clog << "========== cech_complex_from_points ==========" << std::endl;
  Simplex_tree st;
  const int DIMENSION = 3;
  cech_complex_from_points.create_complex(st, DIMENSION);

  // Another way to check num_simplices
  std::clog << "Iterator on Cech complex simplices in the filtration order, with [filtration value]:" << std::endl;
  int num_simplices = 0;
  for (auto f_simplex : st.filtration_simplex_range()) {
    num_simplices++;
    std::clog << "   ( ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> "
              << "[" << st.filtration(f_simplex) << "] ";
    std::clog << std::endl;
  }
  BOOST_CHECK(num_simplices == 15);
  std::clog << "st.num_simplices()=" << st.num_simplices() << std::endl;
  BOOST_CHECK(st.num_simplices() == 15);

  std::clog << "st.dimension()=" << st.dimension() << std::endl;
  BOOST_CHECK(st.dimension() == DIMENSION);
  std::clog << "st.num_vertices()=" << st.num_vertices() << std::endl;
  BOOST_CHECK(st.num_vertices() == 4);

  for (auto f_simplex : st.filtration_simplex_range()) {
    std::clog << "dimension(" << st.dimension(f_simplex) << ") - f = " << st.filtration(f_simplex) << std::endl;
    switch (st.dimension(f_simplex)) {
      case 0:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex), 0.0);
        break;
      case 1:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex), 0.707107, .00001);
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
  double max_radius = 12.0;
  std::clog << "========== OFF FILE NAME = " << off_file_name << " - Cech max_radius=" << max_radius
            << "==========" << std::endl;

  Gudhi::Points_off_reader<Point> off_reader(off_file_name);
  Cech_complex cech_complex_from_file(off_reader.get_point_cloud(), max_radius);

  Simplex_tree stree;
  std::vector<int> simplex = {0, 1, 2};
  stree.insert_simplex_and_subfaces(simplex);
  std::clog << "Check exception throw in debug mode" << std::endl;
  // throw excpt because stree is not empty
  BOOST_CHECK_THROW(cech_complex_from_file.create_complex(stree, 1), std::invalid_argument);
}
#endif
