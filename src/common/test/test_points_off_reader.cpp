/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Points_off_io.h>

#include <iostream>
#include <string>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "points_off_read_write"
#include <boost/test/unit_test.hpp>

using Point_d = std::vector<double>;

BOOST_AUTO_TEST_CASE( points_doc_test )
{
  // Read the OFF file (input file name given as parameter) and triangulates points
  Gudhi::Points_off_reader<Point_d> off_reader("alphacomplexdoc.off");
  // Check the read operation was correct
  BOOST_CHECK(off_reader.is_valid());
  
  // Retrieve the triangulation
  std::vector<Point_d> point_cloud = off_reader.get_point_cloud();
  BOOST_CHECK(point_cloud.size() == 7);

  std::vector<Point_d> expected_points;
  std::vector<double> point = {1.0, 1.0, 0.0};
  expected_points.push_back(Point_d(point.begin(), point.end()));
  point = {7.0, 0.0, 0.0};
  expected_points.push_back(Point_d(point.begin(), point.end()));
  point = {4.0, 6.0, 0.0};
  expected_points.push_back(Point_d(point.begin(), point.end()));
  point = {9.0, 6.0, 0.0};
  expected_points.push_back(Point_d(point.begin(), point.end()));
  point = {0.0, 14.0, 0.0};
  expected_points.push_back(Point_d(point.begin(), point.end()));
  point = {2.0, 19.0, 0.0};
  expected_points.push_back(Point_d(point.begin(), point.end()));
  point = {9.0, 17.0, 0.0};
  expected_points.push_back(Point_d(point.begin(), point.end()));

  BOOST_CHECK(point_cloud == expected_points);
}

BOOST_AUTO_TEST_CASE( Delaunay_triangulation_unexisting_file_read_test )
{
  Gudhi::Points_off_reader<Point_d> off_reader("some_impossible_weird_file_name.off");
  // Check the read operation was correct
  BOOST_CHECK(!off_reader.is_valid());
  
  std::vector<Point_d> point_cloud = off_reader.get_point_cloud();
  BOOST_CHECK(point_cloud.size() == 0);
}
