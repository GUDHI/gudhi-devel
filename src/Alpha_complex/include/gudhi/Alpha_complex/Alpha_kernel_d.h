/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ALPHA_COMPLEX_ALPHA_KERNEL_D_H_
#define ALPHA_COMPLEX_ALPHA_KERNEL_D_H_

#include <CGAL/version.h>  // for CGAL_VERSION_NR

#include <Eigen/Core>  // for EIGEN_VERSION_AT_LEAST

#include <utility>  // for std::make_pair

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if CGAL_VERSION_NR < 1041101000
# error Alpha_complex is only available for CGAL >= 4.11
#endif

#if !EIGEN_VERSION_AT_LEAST(3,1,0)
# error Alpha_complex is only available for Eigen3 >= 3.1.0 installed with CGAL
#endif

namespace Gudhi {

namespace alpha_complex {

/**
 * \class Alpha_kernel_d
 * \brief Alpha complex kernel container.
 * 
 * \details
 * The Alpha complex kernel container stores CGAL Kernel and dispatch basics computations in function of the weighted
 * or not version of the Alpha complex.
 */
template < typename Kernel, bool Weighted = false >
class Alpha_kernel_d {
};

// Unweighted Kernel_d version
template < typename Kernel >
class Alpha_kernel_d<Kernel, false> {
 private:
  // Kernel for functions access.
  Kernel kernel_;
 public:
  // Fake type for compilation to succeed (cf. std::conditional in Alpha_complex.h)
  using Weighted_point_d = void;
  using Point_d = typename Kernel::Point_d;
  // Numeric type of coordinates in the kernel
  using FT = typename Kernel::FT;
  // Sphere is a pair of point and squared radius.
  using Sphere = typename std::pair<Point_d, FT>;

  int get_dimension(const Point_d& p0) const {
    return kernel_.point_dimension_d_object()(p0);
  }

  template<class PointIterator>
  Sphere get_sphere(PointIterator begin, PointIterator end) const {
    Point_d c = kernel_.construct_circumcenter_d_object()(begin, end);
    FT r = kernel_.squared_distance_d_object()(c, *begin);
    return std::make_pair(std::move(c), std::move(r));
  }

  template<class PointIterator>
  FT get_squared_radius(PointIterator begin, PointIterator end) const {
    return kernel_.compute_squared_radius_d_object()(begin, end);
  }

  FT get_squared_radius(const Sphere& sph) const {
    return sph.second;
  }

  bool is_gabriel(const Sphere& circumcenter, const Point_d& point) {
    return kernel_.squared_distance_d_object()(circumcenter.first, point) >= circumcenter.second;
  }
};

// Weighted Kernel_d version
template < typename Kernel >
class Alpha_kernel_d<Kernel, true> {
 private:
  // Kernel for functions access.
  Kernel kernel_;

 public:
  // Fake type for compilation to succeed (cf. std::conditional in Alpha_complex.h)
  using Point_d = void;
  using Weighted_point_d = typename Kernel::Weighted_point_d;
  using Bare_point_d = typename Kernel::Point_d;
  // Numeric type of coordinates in the kernel
  using FT = typename Kernel::FT;
  // Sphere is a weighted point (point + weight [= squared radius]).
  using Sphere = Weighted_point_d;

  int get_dimension(const Weighted_point_d& p0) const {
    return kernel_.point_dimension_d_object()(p0.point());
  }

  template<class PointIterator>
  Sphere get_sphere(PointIterator begin, PointIterator end) const {
    // power_center_d_object has been renamed between CGAL 5.1 and 5.2
#if CGAL_VERSION_NR < 1050200000
    return kernel_.power_center_d_object()(begin, end);
#else
    return kernel_.construct_power_sphere_d_object()(begin, end);
#endif
  }

  template<class PointIterator>
  FT get_squared_radius(PointIterator begin, PointIterator end) const {
    return kernel_.compute_squared_radius_smallest_orthogonal_sphere_d_object()(begin, end);
  }

  FT get_squared_radius(const Sphere& sph) const {
    return sph.weight();
  }

  bool is_gabriel(const Sphere& circumcenter, const Weighted_point_d& point) {
    // power_center_d_object has been renamed between CGAL 5.1 and 5.2
#if CGAL_VERSION_NR < 1050200000
    return kernel_.power_distance_d_object()(circumcenter, point) >= 0;
#else
    return kernel_.compute_power_product_d_object()(circumcenter, point) >= 0;
#endif
  }
};

}  // namespace alpha_complex

namespace alphacomplex = alpha_complex;

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_ALPHA_KERNEL_D_H_