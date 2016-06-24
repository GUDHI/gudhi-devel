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


#include <gudhi/choose_by_farthest_point.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <array>
#include <vector>
#include <iterator>

int main (void)
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4> >   K;
  typedef typename K::FT                           FT;
  typedef typename K::Point_d                      Point_d;
  
  CGAL::Random rd;

  std::vector<Point_d> points;
  for (int i = 0 ; i < 500 ; ++i)
    points.push_back(Point_d(std::array<FT,4>({rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1)})));

  K k;
  std::vector<Point_d> results;
  Gudhi::subsampling::choose_by_farthest_point(k, points, 100, std::back_inserter(results));
  std::cout << "Before sparsification: " << points.size() << " points.\n";
  std::cout << "After  sparsification: " << results.size() << " points.\n";

  return 0;
}
