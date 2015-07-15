/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA Sophia Antipolis-Méditerranée (France)
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

#include <iostream>
#include <vector>

#include <CGAL/Epick_d.h>
#include <CGAL/enum.h>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::Point_d                                 Point_d;
typedef K::FT                                      FT;
typedef K::Hyperplane_d                            Hyperplane_d;
typedef K::Has_on_positive_side_d                  Has_on_positive_side_d;

int main ()
{
  std::vector<Point_d> vertices;
  Point_d v1(std::vector<FT>({-1,1}));
  Point_d v2(std::vector<FT>({1,-1}));
  vertices.push_back(v1);
  vertices.push_back(v2);
  Point_d p(std::vector<FT>({-1,-1}));
  Hyperplane_d hp(vertices.begin(), vertices.end());
  //Hyperplane_d hp(vertices.begin(), vertices.end(), p, CGAL::ON_POSITIVE_SIDE);
  if (Has_on_positive_side_d()(hp, p))
    std::cout << "OK\n";
  else
    std::cout << "NOK\n";
  CGAL::Oriented_side side_p = K::Oriented_side_d()(hp, p);
  if (side_p == CGAL::ZERO)
    std::cout << "Point (-1,-1) is on the line passing through (-1,1) and (1,-1)";
  CGAL::Oriented_side side_v2 = K::Oriented_side_d()(hp, v2);
  if (side_v2 != CGAL::ZERO)
    std::cout << "Point (1,-1) is not on the line passing through (-1,1) and (1,-1)";
}
