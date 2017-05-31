/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA
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

#ifndef DISTANCE_FUNCTIONS_H_
#define DISTANCE_FUNCTIONS_H_

#include <cmath>  // for std::sqrt
#include <type_traits>  // for std::decay
#include <iterator>  // for std::begin, std::end
#include <utility>

namespace Gudhi {

/** @file
 * @brief Global distance functions
 */

/** @brief Compute the Euclidean distance between two Points given by a range of coordinates. The points are assumed to
 * have the same dimension. */
class Euclidean_distance {
 public:
  template< typename Point >
  auto operator()(const Point& p1, const Point& p2) const -> typename std::decay<decltype(*std::begin(p1))>::type {
    auto it1 = p1.begin();
    auto it2 = p2.begin();
    typename Point::value_type dist = 0.;
    for (; it1 != p1.end(); ++it1, ++it2) {
      typename Point::value_type tmp = (*it1) - (*it2);
      dist += tmp*tmp;
    }
    return std::sqrt(dist);
  }
  template< typename T >
  T operator() ( const std::pair< T, T >& f , const std::pair< T, T >& s ) {
    return  sqrt( (f.first-s.first)*(f.first-s.first) + (f.second-s.second)*(f.second-s.second) );
  }
};

}  // namespace Gudhi

#endif  // DISTANCE_FUNCTIONS_H_
