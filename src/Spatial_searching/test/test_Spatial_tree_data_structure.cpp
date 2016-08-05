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
    points.push_back(Point(std::array<FT,4>({rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1)})));

  Points_ds points_ds(points);

  std::size_t closest_pt_index =
    points_ds.query_k_nearest_neighbors(points[10], 1, false).begin()->first;
  BOOST_CHECK(closest_pt_index == 10);

  auto kns_range = points_ds.query_k_nearest_neighbors(points[20], 10, true);

  FT last_dist = -1.;
  for (auto const& nghb : kns_range)
  {
    BOOST_CHECK(nghb.second > last_dist);
    last_dist = nghb.second;
  }

  auto kfs_range = points_ds.query_k_farthest_neighbors(points[20], 10, true);

  last_dist = kfs_range.begin()->second;
  for (auto const& nghb : kfs_range)
  {
    BOOST_CHECK(nghb.second <= last_dist);
    last_dist = nghb.second;
  }
}
