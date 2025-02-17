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
#include <gudhi/MEB_filtration.h>
// to construct Cech_complex from a OFF file of points
#include <gudhi/Points_off_io.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>

#include <CGAL/Epeck_d.h>  // For EXACT or SAFE version

// Type definitions
using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_d;

using Point_cloud = std::vector<Point>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;
using Cech_complex = Gudhi::cech_complex::Cech_complex<Kernel, Simplex_tree>;

BOOST_AUTO_TEST_CASE(Cech_complex_for_documentation) {
  // ----------------------------------------------------------------------------
  //
  // Init of a Cech complex from a point cloud
  //
  // ----------------------------------------------------------------------------
  Point_cloud points;

  std::vector<FT> point0({1., 0.});
  points.emplace_back(point0.begin(), point0.end());
  std::vector<FT> point1({0., 1.});
  points.emplace_back(point1.begin(), point1.end());
  std::vector<FT> point2({2., 1.});
  points.emplace_back(point2.begin(), point2.end());
  std::vector<FT> point3({3., 2.});
  points.emplace_back(point3.begin(), point3.end());
  std::vector<FT> point4({0., 3.});
  points.emplace_back(point4.begin(), point4.end());
  std::vector<FT> point5({3. + std::sqrt(3.), 3.});
  points.emplace_back(point5.begin(), point5.end());
  std::vector<FT> point6({1., 4.});
  points.emplace_back(point6.begin(), point6.end());
  std::vector<FT> point7({3., 4.});
  points.emplace_back(point7.begin(), point7.end());
  std::vector<FT> point8({2., 4. + std::sqrt(3.)});
  points.emplace_back(point8.begin(), point8.end());
  std::vector<FT> point9({0., 4.});
  points.emplace_back(point9.begin(), point9.end());
  std::vector<FT> point10({-0.5, 2.});
  points.emplace_back(point10.begin(), point10.end());

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
      std::clog << ") - distance =" << Gudhi::cech_complex::Sphere_circumradius<Kernel, Filtration_value>()(vp.at(0), vp.at(1))
                << " - filtration =" << st.filtration(f_simplex) << std::endl;
      BOOST_CHECK(vp.size() == 2);
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex),
                                      Gudhi::cech_complex::Sphere_circumradius<Kernel, Filtration_value>()(vp.at(0), vp.at(1)));
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

  Kernel kern;
  Filtration_value f012 = st2.filtration(st2.find({0, 1, 2}));
  std::clog << "f012= " << f012 << std::endl;

  CGAL::NT_converter<FT, Filtration_value> cast_to_fv;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f012, std::sqrt(cast_to_fv(kern.compute_squared_radius_d_object()(points012.begin(), points012.end()))));

  Point_cloud points1410;
  points1410.push_back(cech_complex_for_doc.get_point(1));
  points1410.push_back(cech_complex_for_doc.get_point(4));
  points1410.push_back(cech_complex_for_doc.get_point(10));

  Filtration_value f1410 = st2.filtration(st2.find({1, 4, 10}));
  std::clog << "f1410= " << f1410 << std::endl;

  // In this case, the computed circumsphere using CGAL kernel does not match the minimal enclosing ball; the filtration value check is therefore done against a hardcoded value
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f1410, 1.);

  Point_cloud points469;
  points469.push_back(cech_complex_for_doc.get_point(4));
  points469.push_back(cech_complex_for_doc.get_point(6));
  points469.push_back(cech_complex_for_doc.get_point(9));

  Filtration_value f469 = st2.filtration(st2.find({4, 6, 9}));
  std::clog << "f469= " << f469 << std::endl;

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f469, std::sqrt(cast_to_fv(kern.compute_squared_radius_d_object()(points469.begin(), points469.end()))));

  BOOST_CHECK((st2.find({6, 7, 8}) == st2.null_simplex()));
  BOOST_CHECK((st2.find({3, 5, 7}) == st2.null_simplex()));

  Simplex_tree st3;
  Cech_complex cech5(points, 1e50);
  cech5.create_complex(st3, 5);
  BOOST_CHECK(st3.dimension() == 5);
  auto st3_save = st3;
  std::clog << "Simplex tree from Cech complex\n";
  for (auto f_simplex : st3.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : st3.simplex_vertex_range(f_simplex)) std::clog << vertex << " ";
    std::clog << ") -> " << "[" << st3.filtration(f_simplex) << "]\n";
  }

  st3.reset_filtration(-1); // unnecessary, but ensures we don't cheat
  // Means Output_squared_values=true
  Gudhi::cech_complex::assign_MEB_filtration(Kernel(), st3, points);
  // sqrt all filtration values
  st3.for_each_simplex([&](auto sh, int){
      st3.assign_filtration(sh, std::sqrt(st3.filtration(sh))); });
  std::clog << "Simplex tree from assign_MEB_filtration - after std::sqrt on filtration values\n";
  for (auto f_simplex : st3.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : st3.simplex_vertex_range(f_simplex)) std::clog << vertex << " ";
    std::clog << ") -> " << "[" << st3.filtration(f_simplex) << "]\n";
  }
  BOOST_CHECK(st3 == st3_save); // Should only be an approximate test

  // Same test but with Output_squared_values=false
  st3.reset_filtration(-1); // unnecessary, but ensures we don't cheat
  Gudhi::cech_complex::assign_MEB_filtration<false>(Kernel(), st3, points);
  std::clog << "Simplex tree from assign_MEB_filtration with Output_squared_values=false\n";
  for (auto f_simplex : st3.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : st3.simplex_vertex_range(f_simplex)) std::clog << vertex << " ";
    std::clog << ") -> " << "[" << st3.filtration(f_simplex) << "]\n";
  }
  BOOST_CHECK(st3 == st3_save); // Should only be an approximate test
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
  Filtration_value max_radius = 12.0;
  std::clog << "========== OFF FILE NAME = " << off_file_name << " - Cech max_radius=" << max_radius
            << "==========" << std::endl;

  Gudhi::Points_off_reader<Point> off_reader(off_file_name);
  Cech_complex cech_complex_from_file(off_reader.get_point_cloud(), max_radius);

  Simplex_tree stree;
  std::vector<int> simplex = {0, 1, 2};
  stree.insert_simplex_and_subfaces(simplex);
  std::clog << "Check exception throw in debug mode" << std::endl;
  // throw exception because stree is not empty
  BOOST_CHECK_THROW(cech_complex_from_file.create_complex(stree, 1), std::invalid_argument);
}
#endif
