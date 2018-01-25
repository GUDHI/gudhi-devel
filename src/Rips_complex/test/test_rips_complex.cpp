/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016  INRIA Saclay (France)
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
#define BOOST_TEST_MODULE "rips_complex"
#include <boost/test/unit_test.hpp>

#include <cmath>  // float comparison
#include <limits>
#include <string>
#include <vector>
#include <algorithm>    // std::max

#include <gudhi/Rips_complex.h>
// to construct Rips_complex from a OFF file of points
#include <gudhi/Points_off_io.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Unitary_tests_utils.h>

// Type definitions
using Point = std::vector<double>;
using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Simplex_tree::Filtration_value>;
using Distance_matrix = std::vector<std::vector<Filtration_value>>;

BOOST_AUTO_TEST_CASE(RIPS_DOC_OFF_file) {
  // ----------------------------------------------------------------------------
  //
  // Init of a Rips complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("alphacomplexdoc.off");
  double rips_threshold = 12.0;
  std::cout << "========== OFF FILE NAME = " << off_file_name << " - Rips threshold=" <<
      rips_threshold << "==========" << std::endl;

  Gudhi::Points_off_reader<Point> off_reader(off_file_name);
  Rips_complex rips_complex_from_file(off_reader.get_point_cloud(), rips_threshold, Gudhi::Euclidean_distance());

  const int DIMENSION_1 = 1;
  Simplex_tree st;
  rips_complex_from_file.create_complex(st, DIMENSION_1);
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
  rips_complex_from_file.create_complex(st2, DIMENSION_2);
  std::cout << "st2.dimension()=" << st2.dimension() << std::endl;
  BOOST_CHECK(st2.dimension() == DIMENSION_2);
  
  std::cout << "st2.num_vertices()=" << st2.num_vertices() << std::endl;
  BOOST_CHECK(st2.num_vertices() == NUMBER_OF_VERTICES);

  std::cout << "st2.num_simplices()=" << st2.num_simplices() << std::endl;
  BOOST_CHECK(st2.num_simplices() == 23);

  Simplex_tree::Filtration_value f01 = st2.filtration(st2.find({0, 1}));
  Simplex_tree::Filtration_value f02 = st2.filtration(st2.find({0, 2}));
  Simplex_tree::Filtration_value f12 = st2.filtration(st2.find({1, 2}));
  Simplex_tree::Filtration_value f012 = st2.filtration(st2.find({0, 1, 2}));
  std::cout << "f012= " << f012 << " | f01= " << f01 << " - f02= " << f02 << " - f12= " << f12 << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f012, std::max(f01, std::max(f02,f12)));
  
  Simplex_tree::Filtration_value f45 = st2.filtration(st2.find({4, 5}));
  Simplex_tree::Filtration_value f56 = st2.filtration(st2.find({5, 6}));
  Simplex_tree::Filtration_value f46 = st2.filtration(st2.find({4, 6}));
  Simplex_tree::Filtration_value f456 = st2.filtration(st2.find({4, 5, 6}));
  std::cout << "f456= " << f456 << " | f45= " << f45 << " - f56= " << f56 << " - f46= " << f46 << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f456, std::max(f45, std::max(f56,f46)));

  const int DIMENSION_3 = 3;
  Simplex_tree st3;
  rips_complex_from_file.create_complex(st3, DIMENSION_3);
  std::cout << "st3.dimension()=" << st3.dimension() << std::endl;
  BOOST_CHECK(st3.dimension() == DIMENSION_3);
  
  std::cout << "st3.num_vertices()=" << st3.num_vertices() << std::endl;
  BOOST_CHECK(st3.num_vertices() == NUMBER_OF_VERTICES);

  std::cout << "st3.num_simplices()=" << st3.num_simplices() << std::endl;
  BOOST_CHECK(st3.num_simplices() == 24);

  Simplex_tree::Filtration_value f123 = st3.filtration(st3.find({1, 2, 3}));
  Simplex_tree::Filtration_value f013 = st3.filtration(st3.find({0, 1, 3}));
  Simplex_tree::Filtration_value f023 = st3.filtration(st3.find({0, 2, 3}));
  Simplex_tree::Filtration_value f0123 = st3.filtration(st3.find({0, 1, 2, 3}));
  std::cout << "f0123= " << f0123 << " | f012= " << f012 << " - f123= " << f123 << " - f013= " << f013 <<
      " - f023= " << f023 << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f0123, std::max(f012, std::max(f123, std::max(f013, f023))));

}

using Vector_of_points = std::vector<Point>;

bool is_point_in_list(Vector_of_points points_list, Point point) {
  for (auto& point_in_list : points_list) {
    if (point_in_list == point) {
      return true;  // point found
    }
  }
  return false;  // point not found
}

class Custom_square_euclidean_distance {
 public:
  template< typename Point >
  auto operator()(const Point& p1, const Point& p2) -> typename Point::value_type {
    auto it1 = p1.begin();
    auto it2 = p2.begin();
    typename Point::value_type dist = 0.;
    for (; it1 != p1.end(); ++it1, ++it2) {
      typename Point::value_type tmp = (*it1) - (*it2);
      dist += tmp*tmp;
    }
    return dist;
  }
};

BOOST_AUTO_TEST_CASE(Rips_complex_from_points) {
  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Vector_of_points points;
  std::vector<double> coords = { 0.0, 0.0, 0.0, 1.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 0.0, 0.0, 1.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 0.0, 1.0, 0.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 1.0, 0.0, 0.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));

  // ----------------------------------------------------------------------------
  // Init of a Rips complex from the list of points
  // ----------------------------------------------------------------------------
  Rips_complex rips_complex_from_points(points, 2.0, Custom_square_euclidean_distance());

  std::cout << "========== Rips_complex_from_points ==========" << std::endl;
  Simplex_tree st;
  const int DIMENSION = 3;
  rips_complex_from_points.create_complex(st, DIMENSION);

  // Another way to check num_simplices
  std::cout << "Iterator on Rips complex simplices in the filtration order, with [filtration value]:" << std::endl;
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
      case 2:
      case 3:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex), 2.0);
        break;
      default:
        BOOST_CHECK(false);  // Shall not happen
        break;
    }
  }
}

BOOST_AUTO_TEST_CASE(Rips_doc_csv_file) {
  // ----------------------------------------------------------------------------
  //
  // Init of a Rips complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string csv_file_name("full_square_distance_matrix.csv");
  double rips_threshold = 12.0;
  std::cout << "========== CSV FILE NAME = " << csv_file_name << " - Rips threshold=" <<
      rips_threshold << "==========" << std::endl;

  Distance_matrix distances = Gudhi::read_lower_triangular_matrix_from_csv_file<Filtration_value>(csv_file_name);
  Rips_complex rips_complex_from_file(distances, rips_threshold);

  const int DIMENSION_1 = 1;
  Simplex_tree st;
  rips_complex_from_file.create_complex(st, DIMENSION_1);
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
      std::vector<Simplex_tree::Vertex_handle> vvh;
      std::cout << "vertex = (";
      for (auto vertex : st.simplex_vertex_range(f_simplex)) {
        std::cout << vertex << ",";
        vvh.push_back(vertex);
      }
      std::cout << ") - filtration =" << st.filtration(f_simplex) << std::endl;
      BOOST_CHECK(vvh.size() == 2);
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(f_simplex), distances[vvh.at(0)][vvh.at(1)]);
    }
  }

  const int DIMENSION_2 = 2;
  Simplex_tree st2;
  rips_complex_from_file.create_complex(st2, DIMENSION_2);
  std::cout << "st2.dimension()=" << st2.dimension() << std::endl;
  BOOST_CHECK(st2.dimension() == DIMENSION_2);
  
  std::cout << "st2.num_vertices()=" << st2.num_vertices() << std::endl;
  BOOST_CHECK(st2.num_vertices() == NUMBER_OF_VERTICES);

  std::cout << "st2.num_simplices()=" << st2.num_simplices() << std::endl;
  BOOST_CHECK(st2.num_simplices() == 23);

  Simplex_tree::Filtration_value f01 = st2.filtration(st2.find({0, 1}));
  Simplex_tree::Filtration_value f02 = st2.filtration(st2.find({0, 2}));
  Simplex_tree::Filtration_value f12 = st2.filtration(st2.find({1, 2}));
  Simplex_tree::Filtration_value f012 = st2.filtration(st2.find({0, 1, 2}));
  std::cout << "f012= " << f012 << " | f01= " << f01 << " - f02= " << f02 << " - f12= " << f12 << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f012, std::max(f01, std::max(f02,f12)));
  
  Simplex_tree::Filtration_value f45 = st2.filtration(st2.find({4, 5}));
  Simplex_tree::Filtration_value f56 = st2.filtration(st2.find({5, 6}));
  Simplex_tree::Filtration_value f46 = st2.filtration(st2.find({4, 6}));
  Simplex_tree::Filtration_value f456 = st2.filtration(st2.find({4, 5, 6}));
  std::cout << "f456= " << f456 << " | f45= " << f45 << " - f56= " << f56 << " - f46= " << f46 << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f456, std::max(f45, std::max(f56,f46)));

  const int DIMENSION_3 = 3;
  Simplex_tree st3;
  rips_complex_from_file.create_complex(st3, DIMENSION_3);
  std::cout << "st3.dimension()=" << st3.dimension() << std::endl;
  BOOST_CHECK(st3.dimension() == DIMENSION_3);
  
  std::cout << "st3.num_vertices()=" << st3.num_vertices() << std::endl;
  BOOST_CHECK(st3.num_vertices() == NUMBER_OF_VERTICES);

  std::cout << "st3.num_simplices()=" << st3.num_simplices() << std::endl;
  BOOST_CHECK(st3.num_simplices() == 24);

  Simplex_tree::Filtration_value f123 = st3.filtration(st3.find({1, 2, 3}));
  Simplex_tree::Filtration_value f013 = st3.filtration(st3.find({0, 1, 3}));
  Simplex_tree::Filtration_value f023 = st3.filtration(st3.find({0, 2, 3}));
  Simplex_tree::Filtration_value f0123 = st3.filtration(st3.find({0, 1, 2, 3}));
  std::cout << "f0123= " << f0123 << " | f012= " << f012 << " - f123= " << f123 << " - f013= " << f013 <<
      " - f023= " << f023 << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(f0123, std::max(f012, std::max(f123, std::max(f013, f023))));

}

#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE(Rips_create_complex_throw) {
  // ----------------------------------------------------------------------------
  //
  // Init of a Rips complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("alphacomplexdoc.off");
  double rips_threshold = 12.0;
  std::cout << "========== OFF FILE NAME = " << off_file_name << " - Rips threshold=" <<
      rips_threshold << "==========" << std::endl;

  Gudhi::Points_off_reader<Point> off_reader(off_file_name);
  Rips_complex rips_complex_from_file(off_reader.get_point_cloud(), rips_threshold, Gudhi::Euclidean_distance());

  Simplex_tree stree;
  std::vector<int> simplex = {0, 1, 2};
  stree.insert_simplex_and_subfaces(simplex);
  std::cout << "Check exception throw in debug mode" << std::endl;
  // throw excpt because stree is not empty
  BOOST_CHECK_THROW (rips_complex_from_file.create_complex(stree, 1), std::invalid_argument);
}
#endif
