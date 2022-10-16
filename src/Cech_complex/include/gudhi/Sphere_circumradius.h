/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hind Montassif
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SPHERE_CIRCUMRADIUS_H_
#define SPHERE_CIRCUMRADIUS_H_

#include <CGAL/Epick_d.h> // for #include <CGAL/NT_converter.h> which is not working/compiling alone
#include <CGAL/Lazy_exact_nt.h> // for CGAL::exact

#include <cmath>  // for std::sqrt
#include <vector>

namespace Gudhi {

namespace cech_complex {

/** \private @brief Compute the circumradius of the sphere passing through points given by a range of coordinates.
 * The points are assumed to have the same dimension. */
template<typename Kernel, typename Filtration_value>
class Sphere_circumradius {
 private:
    Kernel kernel_;
    const bool exact_;
 public:
    using FT = typename Kernel::FT;
    using Point = typename Kernel::Point_d;
    using Point_cloud = typename std::vector<Point>;

    CGAL::NT_converter<FT, Filtration_value> cast_to_fv;

   /** \brief Circumradius of sphere passing through two points using CGAL.
   *
   * @param[in] point_1
   * @param[in] point_2
   * @return Sphere circumradius passing through two points.
   * \tparam Point must be a Kernel::Point_d from CGAL.
   *
   */
  Filtration_value operator()(const Point& point_1, const Point& point_2) const {
    auto squared_dist_obj = kernel_.squared_distance_d_object()(point_1, point_2);
    if(exact_) CGAL::exact(squared_dist_obj);
    return std::sqrt(cast_to_fv(squared_dist_obj)) / 2.;
  }

  /** \brief Circumradius of sphere passing through point cloud using CGAL.
   *
   * @param[in] point_cloud The points.
   * @return Sphere circumradius passing through the points.
   * \tparam Point_cloud must be a range of Kernel::Point_d points from CGAL.
   *
   */
  Filtration_value operator()(const Point_cloud& point_cloud) const {
    auto squared_radius_obj = kernel_.compute_squared_radius_d_object()(point_cloud.begin(), point_cloud.end());
    if(exact_) CGAL::exact(squared_radius_obj);
    return std::sqrt(cast_to_fv(squared_radius_obj));
  }

  /** \brief Constructor
   * @param[in] exact Option for exact filtration values computation. Not exact if `Kernel` is not <a target="_blank"
   * href="https://doc.cgal.org/latest/Kernel_d/structCGAL_1_1Epeck__d.html">CGAL::Epeck_d</a>.
   * Default is false.
   */
  Sphere_circumradius(const bool exact = false) : exact_(exact) {}

};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // SPHERE_CIRCUMRADIUS_H_
