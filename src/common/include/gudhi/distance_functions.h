/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DISTANCE_FUNCTIONS_H_
#define DISTANCE_FUNCTIONS_H_

#include <gudhi/Debug_utils.h>

#include <gudhi/Miniball.hpp>

#include <boost/range/metafunctions.hpp>
#include <boost/range/size.hpp>

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

/** @brief Compute the radius of the minimal enclosing ball between Points given by a range of coordinates.
 * The points are assumed to have the same dimension. */
class Minimal_enclosing_ball_radius {
 public:
  /** \brief Minimal_enclosing_ball_radius from two points.
   *
   * @param[in] point_1 First point.
   * @param[in] point_2 second point.
   * @return The minimal enclosing ball radius for the two points (aka. Euclidean distance / 2.).
   *
   * \tparam Point must be a range of Cartesian coordinates.
   *
   */
  template< typename Point >
  typename std::iterator_traits<typename boost::range_iterator<Point>::type>::value_type
  operator()(const Point& point_1, const Point& point_2) const {
    return Euclidean_distance()(point_1, point_2) / 2.;
  }
  /** \brief Minimal_enclosing_ball_radius from a point cloud.
   *
   * @param[in] point_cloud The points.
   * @return The minimal enclosing ball radius for the points.
   *
   * \tparam Point_cloud must be a range of points with Cartesian coordinates.
   * Point_cloud is a range over a range of Coordinate.
   *
   */
  template< typename Point_cloud,
            typename Point_iterator = typename boost::range_const_iterator<Point_cloud>::type,
            typename Point = typename std::iterator_traits<Point_iterator>::value_type,
            typename Coordinate_iterator = typename boost::range_const_iterator<Point>::type,
            typename Coordinate = typename std::iterator_traits<Coordinate_iterator>::value_type>
  Coordinate
  operator()(const Point_cloud& point_cloud) const {
    using Min_sphere = Miniball::Miniball<Miniball::CoordAccessor<Point_iterator, Coordinate_iterator>>;

    Min_sphere ms(boost::size(*point_cloud.begin()), point_cloud.begin(), point_cloud.end());
#ifdef DEBUG_TRACES
    std::cout << "Minimal_enclosing_ball_radius = " << std::sqrt(ms.squared_radius()) << " | nb points = "
              << boost::size(point_cloud) << " | dimension = "
              << boost::size(*point_cloud.begin()) << std::endl;
#endif  // DEBUG_TRACES

    return std::sqrt(ms.squared_radius());
  }
};

}  // namespace Gudhi

#endif  // DISTANCE_FUNCTIONS_H_
