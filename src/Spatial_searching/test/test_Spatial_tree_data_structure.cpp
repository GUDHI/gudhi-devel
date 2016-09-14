/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016  INRIA Sophia-Antipolis (France)
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
#define BOOST_TEST_MODULE Spatial_searching - test Spatial_tree_data_structure
#include <boost/test/unit_test.hpp>

#include <gudhi/Spatial_tree_data_structure.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <array>
#include <vector>

BOOST_AUTO_TEST_CASE(test_Spatial_tree_data_structure)
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4> >    K;
  typedef K::FT                                     FT;
  typedef K::Point_d                                Point;
  typedef std::vector<Point>                        Points;

  typedef Gudhi::spatial_searching::Spatial_tree_data_structure<
    K, Points>                                      Points_ds;
  
  CGAL::Random rd;

  Points points;
  for (int i = 0 ; i < 500 ; ++i)
    points.push_back(Point(rd.get_double(-1.,1), rd.get_double(-1.,1), rd.get_double(-1.,1), rd.get_double(-1.,1)));

  Points_ds points_ds(points);

  // Test query_k_nearest_neighbors
  std::size_t closest_pt_index =
    points_ds.query_k_nearest_neighbors(points[10], 1, false).begin()->first;
  BOOST_CHECK(closest_pt_index == 10);

  auto kns_range = points_ds.query_k_nearest_neighbors(points[20], 10, true);

  std::vector<std::size_t> knn_result;
  FT last_dist = -1.;
  for (auto const& nghb : kns_range)
  {
    BOOST_CHECK(nghb.second > last_dist);
    knn_result.push_back(nghb.second);
    last_dist = nghb.second;
  }

  // Test query_incremental_nearest_neighbors 
  closest_pt_index =
    points_ds.query_incremental_nearest_neighbors(points[10]).begin()->first;
  BOOST_CHECK(closest_pt_index == 10);

  auto ins_range = points_ds.query_incremental_nearest_neighbors(points[20]);

  std::vector<std::size_t> inn_result;
  last_dist = -1.;
  auto ins_it = ins_range.begin();
  for (int i = 0 ; i < 10 ; ++ins_it, ++i)
  {
    auto const& nghb = *ins_it;
    BOOST_CHECK(nghb.second > last_dist);
    inn_result.push_back(nghb.second);
    last_dist = nghb.second;
  }

  // Same result for KNN and INN?
  BOOST_CHECK(knn_result == inn_result);

  // Test query_k_farthest_neighbors
  auto kfs_range = points_ds.query_k_farthest_neighbors(points[20], 10, true);

  std::vector<std::size_t> kfn_result;
  last_dist = kfs_range.begin()->second;
  for (auto const& nghb : kfs_range)
  {
    BOOST_CHECK(nghb.second <= last_dist);
    kfn_result.push_back(nghb.second);
    last_dist = nghb.second;
  }

  // Test query_k_farthest_neighbors
  auto ifs_range = points_ds.query_incremental_farthest_neighbors(points[20]);

  std::vector<std::size_t> ifn_result;
  last_dist = ifs_range.begin()->second;
  auto ifs_it = ifs_range.begin();
  for (int i = 0; i < 10; ++ifs_it, ++i)
  {
    auto const& nghb = *ifs_it;
    BOOST_CHECK(nghb.second <= last_dist);
    ifn_result.push_back(nghb.second);
    last_dist = nghb.second;
  }

  // Same result for KFN and IFN?
  BOOST_CHECK(kfn_result == ifn_result);
}
