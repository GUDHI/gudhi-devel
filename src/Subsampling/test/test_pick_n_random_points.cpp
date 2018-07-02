/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016 Inria
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
#define BOOST_TEST_MODULE Subsampling - test pick_n_random_points
#include <boost/test/unit_test.hpp>

#include <gudhi/pick_n_random_points.h>
#include <vector>
#include <iterator>

#include <CGAL/Epick_d.h>


BOOST_AUTO_TEST_CASE(test_pick_n_random_points)
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>                K;
  typedef typename K::FT                                            FT;
  typedef typename K::Point_d                                       Point_d;
  
  std::vector<Point_d> vect;
  vect.push_back(Point_d(std::vector<FT>({0,0,0,0})));
  vect.push_back(Point_d(std::vector<FT>({0,0,0,1})));
  vect.push_back(Point_d(std::vector<FT>({0,0,1,0})));
  vect.push_back(Point_d(std::vector<FT>({0,0,1,1})));
  vect.push_back(Point_d(std::vector<FT>({0,1,0,0})));
  vect.push_back(Point_d(std::vector<FT>({0,1,0,1})));
  vect.push_back(Point_d(std::vector<FT>({0,1,1,0})));
  vect.push_back(Point_d(std::vector<FT>({0,1,1,1})));
  vect.push_back(Point_d(std::vector<FT>({1,0,0,0})));
  vect.push_back(Point_d(std::vector<FT>({1,0,0,1})));
  vect.push_back(Point_d(std::vector<FT>({1,0,1,0})));
  vect.push_back(Point_d(std::vector<FT>({1,0,1,1})));
  vect.push_back(Point_d(std::vector<FT>({1,1,0,0})));
  vect.push_back(Point_d(std::vector<FT>({1,1,0,1})));
  vect.push_back(Point_d(std::vector<FT>({1,1,1,0})));
  vect.push_back(Point_d(std::vector<FT>({1,1,1,1})));

  std::vector<Point_d> results;
  Gudhi::subsampling::pick_n_random_points(vect, 5, std::back_inserter(results));
  std::cout << "landmark vector contains: ";
  for (auto l: results)
    std::cout << l << "\n";
  
  BOOST_CHECK(results.size() == 5);
}
