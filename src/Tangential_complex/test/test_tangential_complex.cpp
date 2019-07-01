/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Tangential_complex - test tangential complex
#include <boost/test/unit_test.hpp>

#include <gudhi/Tangential_complex.h>
#include <gudhi/sparsify_point_set.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <array>
#include <vector>

namespace tc = Gudhi::tangential_complex;

BOOST_AUTO_TEST_CASE(test_Spatial_tree_data_structure) {
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
  typedef Kernel::Point_d Point;
  typedef tc::Tangential_complex<
      Kernel, CGAL::Dynamic_dimension_tag,
      CGAL::Parallel_tag> TC;

  const int INTRINSIC_DIM = 2;
  const int AMBIENT_DIM = 3;
  const int NUM_POINTS = 50;

  Kernel k;

  // Generate points on a 2-sphere
  CGAL::Random_points_on_sphere_d<Point> generator(AMBIENT_DIM, 3.);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0; i < NUM_POINTS; ++i)
    points.push_back(*generator++);

  // Compute the TC
  TC tc(points, INTRINSIC_DIM, k);
  tc.compute_tangential_complex();

  // Try to fix inconsistencies. Give it 60 seconds to succeed
  auto perturb_ret = tc.fix_inconsistencies_using_perturbation(0.01, 60);

  BOOST_CHECK(perturb_ret.success);

  // Export the TC into a Simplex_tree
  Gudhi::Simplex_tree<> stree;
  tc.create_complex(stree);
}

BOOST_AUTO_TEST_CASE(test_mini_tangential) {
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
  typedef Kernel::Point_d Point;
  typedef tc::Tangential_complex<Kernel, CGAL::Dynamic_dimension_tag, CGAL::Parallel_tag> TC;


  const int INTRINSIC_DIM = 1;

  // Generate points on a 2-sphere
  std::vector<Point> points;
  // [[0, 0], [1, 0], [0, 1], [1, 1]]
  std::vector<double> point = {0.0, 0.0};
  points.push_back(Point(point.size(), point.begin(), point.end()));
  point = {1.0, 0.0};
  points.push_back(Point(point.size(), point.begin(), point.end()));
  point = {0.0, 1.0};
  points.push_back(Point(point.size(), point.begin(), point.end()));
  point = {1.0, 1.0};
  points.push_back(Point(point.size(), point.begin(), point.end()));
  std::cout << "points = " << points.size() << std::endl;
  Kernel k;

  // Compute the TC
  TC tc(points, INTRINSIC_DIM, k);
  tc.compute_tangential_complex();
  TC::Num_inconsistencies num_inc = tc.number_of_inconsistent_simplices();
  std::cout << "TC vertices = " << tc.number_of_vertices() << " - simplices = " << num_inc.num_simplices <<
               " - inc simplices = " << num_inc.num_inconsistent_simplices <<
               " - inc stars = " << num_inc.num_inconsistent_stars << std::endl;

  BOOST_CHECK(tc.number_of_vertices() == 4);
  BOOST_CHECK(num_inc.num_simplices == 4);
  BOOST_CHECK(num_inc.num_inconsistent_simplices == 0);
  BOOST_CHECK(num_inc.num_inconsistent_stars == 0);

  // Export the TC into a Simplex_tree
  Gudhi::Simplex_tree<> stree;
  tc.create_complex(stree);
  std::cout << "ST vertices = " << stree.num_vertices() << " - simplices = " << stree.num_simplices() << std::endl;

  BOOST_CHECK(stree.num_vertices() == 4);
  BOOST_CHECK(stree.num_simplices() == 6);

  tc.fix_inconsistencies_using_perturbation(0.01, 30.0);

  BOOST_CHECK(tc.number_of_vertices() == 4);
  BOOST_CHECK(num_inc.num_simplices == 4);
  BOOST_CHECK(num_inc.num_inconsistent_simplices == 0);
  BOOST_CHECK(num_inc.num_inconsistent_stars == 0);

  // Export the TC into a Simplex_tree
  tc.create_complex(stree);
  std::cout << "ST vertices = " << stree.num_vertices() << " - simplices = " << stree.num_simplices() << std::endl;

  BOOST_CHECK(stree.num_vertices() == 4);
  BOOST_CHECK(stree.num_simplices() == 6);
}

#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE(test_basic_example_throw) {
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
  typedef Kernel::FT FT;
  typedef Kernel::Point_d Point;
  typedef Kernel::Vector_d Vector;
  typedef tc::Tangential_complex<Kernel, CGAL::Dynamic_dimension_tag,CGAL::Parallel_tag> TC;

  const int INTRINSIC_DIM = 2;
  const int AMBIENT_DIM = 3;
  const int NUM_POINTS = 1000;

  Kernel k;

  // Generate points on a 2-sphere
  CGAL::Random_points_on_sphere_d<Point> generator(AMBIENT_DIM, 3.);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0; i < NUM_POINTS; ++i)
    points.push_back(*generator++);

  // Compute the TC
  TC tc(points, INTRINSIC_DIM, k);
  tc.set_max_squared_edge_length(0.01);
  std::cout << "test_basic_example_throw - set_max_squared_edge_length(0.01) to make GUDHI_CHECK fail" << std::endl;
  BOOST_CHECK_THROW(tc.compute_tangential_complex(), std::invalid_argument);

}
#endif
