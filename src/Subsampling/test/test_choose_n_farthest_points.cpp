/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

// #ifdef _DEBUG
// # define TBB_USE_THREADING_TOOL
// #endif

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Subsampling - test choose_n_farthest_points
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/choose_n_farthest_points.h>
#include <vector>
#include <iterator>

#include <CGAL/Epick_d.h>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef typename K::FT FT;
typedef typename K::Point_d Point_d;

typedef boost::mpl::list<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>, CGAL::Epick_d<CGAL::Dimension_tag<4>>> list_of_tested_kernels;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_choose_farthest_point, Kernel, list_of_tested_kernels) {
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_d Point_d;
  std::vector< Point_d > points, landmarks;
  // Add grid points (625 points)
  for (FT i = 0; i < 5; i += 1.0)
    for (FT j = 0; j < 5; j += 1.0)
      for (FT k = 0; k < 5; k += 1.0)
        for (FT l = 0; l < 5; l += 1.0) {
          std::vector<FT> point({i, j, k, l});
          points.emplace_back(point.begin(), point.end());
        }

  landmarks.clear();
  Kernel k;
  Gudhi::subsampling::choose_n_farthest_points(k, points, 100, Gudhi::subsampling::random_starting_point, std::back_inserter(landmarks));

  BOOST_CHECK(landmarks.size() == 100);
  for (auto landmark : landmarks)
  {
    // Check all landmarks are in points
    BOOST_CHECK(std::find (points.begin(), points.end(), landmark) != points.end());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_choose_farthest_point_limits, Kernel, list_of_tested_kernels) {
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_d Point_d;
  std::vector< Point_d > points, landmarks;
  std::vector< FT > distances;
  landmarks.clear();
  Kernel k;
  // Choose -1 farthest points in an empty point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, -1, -1, std::back_inserter(landmarks), std::back_inserter(distances));
  BOOST_CHECK(landmarks.size() == 0);
  landmarks.clear(); distances.clear();
  // Choose 0 farthest points in an empty point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, 0, -1, std::back_inserter(landmarks), std::back_inserter(distances));
  BOOST_CHECK(landmarks.size() == 0);
  landmarks.clear(); distances.clear();
  // Choose 1 farthest points in an empty point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, 1, -1, std::back_inserter(landmarks), std::back_inserter(distances));
  BOOST_CHECK(landmarks.size() == 0);
  landmarks.clear(); distances.clear();

  std::vector<FT> point({0.0, 0.0, 0.0, 0.0});
  points.emplace_back(point.begin(), point.end());
  // Choose -1 farthest points in a one point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, -1, -1, std::back_inserter(landmarks), std::back_inserter(distances));
  BOOST_CHECK(landmarks.size() == 1 && distances.size() == 1);
  BOOST_CHECK(distances[0] == std::numeric_limits<FT>::infinity());
  landmarks.clear(); distances.clear();
  // Choose 0 farthest points in a one point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, 0, -1, std::back_inserter(landmarks), std::back_inserter(distances));
  BOOST_CHECK(landmarks.size() == 0 && distances.size() == 0);
  landmarks.clear(); distances.clear();
  // Choose 1 farthest points in a one point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, 1, -1, std::back_inserter(landmarks), std::back_inserter(distances));
  BOOST_CHECK(landmarks.size() == 1 && distances.size() == 1);
  BOOST_CHECK(distances[0] == std::numeric_limits<FT>::infinity());
  landmarks.clear(); distances.clear();

  std::vector<FT> point2({1.0, 0.0, 0.0, 0.0});
  points.emplace_back(point2.begin(), point2.end());
  // Choose all farthest points among 2 points
  Gudhi::subsampling::choose_n_farthest_points(k, points, -1, -1, std::back_inserter(landmarks), std::back_inserter(distances));
  BOOST_CHECK(landmarks.size() == 2 && distances.size() == 2);
  BOOST_CHECK(distances[0] == std::numeric_limits<FT>::infinity());
  BOOST_CHECK(distances[1] == 1);
  landmarks.clear(); distances.clear();

  // Ignore duplicated points
  points.emplace_back(point.begin(), point.end());
  Gudhi::subsampling::choose_n_farthest_points(k, points, -1, -1, std::back_inserter(landmarks), std::back_inserter(distances));
  BOOST_CHECK(landmarks.size() == 2 && distances.size() == 2);
  BOOST_CHECK(distances[0] == std::numeric_limits<FT>::infinity());
  BOOST_CHECK(distances[1] == 1);
  landmarks.clear(); distances.clear();
}
