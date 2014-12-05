 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       Clément Maria
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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

/* Compute the Euclidean distance between two Points given
 * by a range of coordinates. The points are assumed to have 
 * the same dimension. */
template< typename Point >
double euclidean_distance( Point &p1, Point &p2)
{
  double dist = 0.;
  auto it1 = p1.begin(); auto it2 = p2.begin();
  for(; it1 != p1.end(); ++it1, ++it2) 
    { 
      double tmp = *it1 - *it2;
      dist += tmp*tmp; 
    }
    return sqrt(dist);
}
