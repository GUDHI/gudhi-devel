/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016  INRIA
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
#define BOOST_TEST_MODULE Spatial_searching - test Kd_tree_search
#include <boost/test/unit_test.hpp>

#include <gudhi/Kd_tree_search.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <vector>

BOOST_AUTO_TEST_CASE(test_Kd_tree_search) {
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4> > K;
  typedef K::FT FT;
  typedef K::Point_d Point;
  typedef std::vector<Point> Points;

  typedef Gudhi::spatial_searching::Kd_tree_search<
      K, Points> Points_ds;

  CGAL::Random rd;

  Points points;
  for (int i = 0; i < 500; ++i)
    points.push_back(Point(rd.get_double(-1., 1), rd.get_double(-1., 1), rd.get_double(-1., 1), rd.get_double(-1., 1)));

  Points_ds points_ds(points);

  // Test k_nearest_neighbors
  std::size_t closest_pt_index =
      points_ds.k_nearest_neighbors(points[10], 1, false).begin()->first;
  BOOST_CHECK(closest_pt_index == 10);

  auto kns_range = points_ds.k_nearest_neighbors(points[20], 10, true);

  std::vector<std::size_t> knn_result;
  FT last_dist = -1.;
  for (auto const& nghb : kns_range) {
    BOOST_CHECK(nghb.second > last_dist);
    knn_result.push_back(nghb.second);
    last_dist = nghb.second;
  }

  // Test incremental_nearest_neighbors 
  closest_pt_index =
      points_ds.incremental_nearest_neighbors(points[10]).begin()->first;
  BOOST_CHECK(closest_pt_index == 10);

  auto inn_range = points_ds.incremental_nearest_neighbors(points[20]);

  std::vector<std::size_t> inn_result;
  last_dist = -1.;
  auto inn_it = inn_range.begin();
  for (int i = 0; i < 10; ++inn_it, ++i) {
    auto const& nghb = *inn_it;
    BOOST_CHECK(nghb.second > last_dist);
    inn_result.push_back(nghb.second);
    last_dist = nghb.second;
  }

  // Same result for KNN and INN?
  BOOST_CHECK(knn_result == inn_result);

  // Test k_furthest_neighbors
  auto kfn_range = points_ds.k_furthest_neighbors(points[20], 10, true);

  std::vector<std::size_t> kfn_result;
  last_dist = kfn_range.begin()->second;
  for (auto const& nghb : kfn_range) {
    BOOST_CHECK(nghb.second <= last_dist);
    kfn_result.push_back(nghb.second);
    last_dist = nghb.second;
  }

  // Test k_furthest_neighbors
  auto ifn_range = points_ds.incremental_furthest_neighbors(points[20]);

  std::vector<std::size_t> ifn_result;
  last_dist = ifn_range.begin()->second;
  auto ifn_it = ifn_range.begin();
  for (int i = 0; i < 10; ++ifn_it, ++i) {
    auto const& nghb = *ifn_it;
    BOOST_CHECK(nghb.second <= last_dist);
    ifn_result.push_back(nghb.second);
    last_dist = nghb.second;
  }

  // Same result for KFN and IFN?
  BOOST_CHECK(kfn_result == ifn_result);

  // Test all_near_neighbors
  Point rs_q(rd.get_double(-1., 1), rd.get_double(-1., 1), rd.get_double(-1., 1), rd.get_double(-1., 1));
  std::vector<std::size_t> rs_result;
  points_ds.all_near_neighbors(rs_q, 0.5, std::back_inserter(rs_result));
  K k;
  for (auto const& p_idx : rs_result)
    BOOST_CHECK(k.squared_distance_d_object()(points[p_idx], rs_q) <= 0.5);
}
