/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hind Montassif
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CECH_KERNEL_H_
#define CECH_KERNEL_H_

#include <CGAL/Epeck_d.h>

#include <cmath>  // for std::sqrt

namespace Gudhi {

// namespace cech_complex {

/** @brief Compute the radius of the minimal enclosing ball between Points given by a range of coordinates.
 * The points are assumed to have the same dimension. */
class Minimal_enclosing_ball_radius {
 public:
   /** \brief Enclosing ball radius from two points using CGAL.
   *
   * @param[in] point_1
   * @param[in] point_2
   * @return Enclosing ball radius for the two points.
   * \tparam Point must be a Kernel::Point_d from CGAL.
   *
   */
  template< typename Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>,
            typename Point= typename Kernel::Point_d>
  double operator()(const Point& point_1, const Point& point_2) const {
    Kernel kernel_;
    return std::sqrt(CGAL::to_double(kernel_.squared_distance_d_object()(point_1, point_2))) / 2.;
  }

  /** \brief Minimal_enclosing_ball_radius from two points.
   *
   * @param[in] point_1 First point.
   * @param[in] point_2 second point.
   * @return The minimal enclosing ball radius for the two points (aka. Euclidean distance / 2.).
   *
   * \tparam Point must be a range of Cartesian coordinates.
   *
   */
//   template< typename Point >
//   typename std::iterator_traits<typename boost::range_iterator<Point>::type>::value_type
//   operator()(const Point& point_1, const Point& point_2) const {
//     std::clog << "#" << *point_1.begin() << "##" << *point_2.begin() << std::endl;
//     return Euclidean_distance()(point_1, point_2) / 2.;
//   }

  /** \brief Enclosing ball radius from a point cloud using CGAL.
   *
   * @param[in] point_cloud The points.
   * @return Enclosing ball radius for the points.
   * \tparam Point_cloud must be a range of Kernel::Point_d points from CGAL.
   *
   */
  template< typename Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>,
            typename Point= typename Kernel::Point_d,
            typename Point_cloud = std::vector<Point>>
  double operator()(const Point_cloud& point_cloud) const {
    Kernel kernel_;
    return std::sqrt(CGAL::to_double(kernel_.compute_squared_radius_d_object()(point_cloud.begin(), point_cloud.end())));
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
//   template< typename Point_cloud,
//             typename Point_iterator = typename boost::range_const_iterator<Point_cloud>::type,
//             typename Point = typename std::iterator_traits<Point_iterator>::value_type,
//             typename Coordinate_iterator = typename boost::range_const_iterator<Point>::type,
//             typename Coordinate = typename std::iterator_traits<Coordinate_iterator>::value_type>
//   Coordinate
//   operator()(const Point_cloud& point_cloud) const {
//     using Min_sphere = Miniball::Miniball<Miniball::CoordAccessor<Point_iterator, Coordinate_iterator>>;
// 
//     Min_sphere ms(boost::size(*point_cloud.begin()), point_cloud.begin(), point_cloud.end());
// #ifdef DEBUG_TRACES
//     std::clog << "Minimal_enclosing_ball_radius = " << std::sqrt(ms.squared_radius()) << " | nb points = "
//               << boost::size(point_cloud) << " | dimension = "
//               << boost::size(*point_cloud.begin()) << std::endl;
// #endif  // DEBUG_TRACES
// 
//     return std::sqrt(ms.squared_radius());
//   }
};

/**
 * \class Cech_kernel
 * \brief Cech complex kernel container.
 * 
 * \details
 * The Cech complex kernel container stores CGAL Kernel and dispatch basic computations.
 */

// template < typename Kernel >
// class Cech_kernel<Kernel> {
//  private:
//   // Kernel for functions access.
//   Kernel kernel_;
//  public:
//   using Point_d = typename Kernel::Point_d;
//   // Numeric type of coordinates in the kernel
//   using FT = typename Kernel::FT;
//   // Sphere is a pair of point and squared radius.
//   using Sphere = typename std::pair<Point_d, FT>;
// 
//   int get_dimension(const Point_d& p0) const {
//     return kernel_.point_dimension_d_object()(p0);
//   }
// 
//   template<class PointIterator>
//   Sphere get_sphere(PointIterator begin, PointIterator end) const {
//     Point_d c = kernel_.construct_circumcenter_d_object()(begin, end);
//     FT r = kernel_.squared_distance_d_object()(c, *begin);
//     return std::make_pair(std::move(c), std::move(r));
//   }
// 
//   template<class PointIterator>
//   FT get_squared_radius(PointIterator begin, PointIterator end) const {
//     return kernel_.compute_squared_radius_d_object()(begin, end);
//   }
// 
//   FT get_squared_radius(const Sphere& sph) const {
//     return sph.second;
//   }
// };


//}  // namespace cech_complex

// namespace cechcomplex = cech_complex;

}  // namespace Gudhi

#endif  // CECH_KERNEL_H_