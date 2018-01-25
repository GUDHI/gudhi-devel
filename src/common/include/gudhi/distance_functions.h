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

#include <gudhi/Debug_utils.h>

#include <boost/range/metafunctions.hpp>

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
  // boost::range_value is not SFINAE-friendly so we cannot use it in the return type
  template< typename Point >
  typename std::iterator_traits<typename boost::range_iterator<Point>::type>::value_type
  operator()(const Point& p1, const Point& p2) const {
    auto it1 = std::begin(p1);
    auto it2 = std::begin(p2);
    typedef typename boost::range_value<Point>::type NT;
    NT dist = 0;
    for (; it1 != std::end(p1); ++it1, ++it2) {
      GUDHI_CHECK(it2 != std::end(p2), "inconsistent point dimensions");
      NT tmp = *it1 - *it2;
      dist += tmp*tmp;
    }
    GUDHI_CHECK(it2 == std::end(p2), "inconsistent point dimensions");
    using std::sqrt;
    return sqrt(dist);
  }
  template< typename T >
  T operator() (const std::pair< T, T >& f, const std::pair< T, T >& s) const {
    T dx = f.first - s.first;
    T dy = f.second - s.second;
    using std::sqrt;
    return sqrt(dx*dx + dy*dy);
  }
};

}  // namespace Gudhi

#endif  // DISTANCE_FUNCTIONS_H_
