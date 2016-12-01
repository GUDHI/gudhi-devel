/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016 INRIA
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

// #ifdef _DEBUG
// # define TBB_USE_THREADING_TOOL
// #endif

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "witness_complex_points"
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
          points.push_back(Point_d(point.begin(), point.end()));
        }

  landmarks.clear();
  Kernel k;
  Gudhi::subsampling::choose_n_farthest_points(k, points, 100, std::back_inserter(landmarks));

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
  landmarks.clear();
  Kernel k;
  // Choose -1 farthest points in an empty point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, -1, std::back_inserter(landmarks));
  BOOST_CHECK(landmarks.size() == 0);
  landmarks.clear();
  // Choose 0 farthest points in an empty point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, 0, std::back_inserter(landmarks));
  BOOST_CHECK(landmarks.size() == 0);
  landmarks.clear();
  // Choose 1 farthest points in an empty point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, 1, std::back_inserter(landmarks));
  BOOST_CHECK(landmarks.size() == 0);
  landmarks.clear();

  std::vector<FT> point({0.0, 0.0, 0.0, 0.0});
  points.push_back(Point_d(point.begin(), point.end()));
  // Choose -1 farthest points in an empty point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, -1, std::back_inserter(landmarks));
  BOOST_CHECK(landmarks.size() == 1);
  landmarks.clear();
  // Choose 0 farthest points in a one point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, 0, std::back_inserter(landmarks));
  BOOST_CHECK(landmarks.size() == 0);
  landmarks.clear();
  // Choose 1 farthest points in a one point cloud
  Gudhi::subsampling::choose_n_farthest_points(k, points, 1, std::back_inserter(landmarks));
  BOOST_CHECK(landmarks.size() == 1);
  landmarks.clear();

}
