/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <boost/test/unit_test.hpp>

#include <CGAL/NT_converter.h>

#include <iostream>
#include <vector>
#include <utility>  // for std::pair

#include <gudhi/Alpha_complex/Alpha_kernel_d.h>
#include <gudhi/Unitary_tests_utils.h>


template<class CGAL_kernel>
void test_alpha_kernel_d_dimension() {
  // Test for a point (weighted or not) in 4d, that the dimension is 4.

  Gudhi::alpha_complex::Alpha_kernel_d<CGAL_kernel, false> kernel;
  std::vector<double> p0 {0., 1., 2., 3.};
  typename CGAL_kernel::Point_d p0_d(p0.begin(), p0.end());

  std::clog << "Dimension is " << kernel.get_dimension(p0_d) << std::endl;
  BOOST_CHECK(kernel.get_dimension(p0_d) == 4);

  Gudhi::alpha_complex::Alpha_kernel_d<CGAL_kernel, true> w_kernel;
  typename CGAL_kernel::Weighted_point_d w_p0_d(p0_d, 10.);

  std::clog << "Dimension is " << w_kernel.get_dimension(w_p0_d) << std::endl;
  BOOST_CHECK(w_kernel.get_dimension(w_p0_d) == 4);
}

template<class CGAL_kernel>
void test_alpha_kernel_d_get_sphere() {
  // Test with 5 points on a 3-sphere, that get_sphere returns the same center and squared radius
  // for dD unweighted and for dD weighted with all weights at 0.

  using Unweighted_kernel = Gudhi::alpha_complex::Alpha_kernel_d<CGAL_kernel, false>;
  // Sphere: (x-1)² + (y-1)² + z² + t² = 1
  // At least 5 points for a 3-sphere
  std::vector<double> p0 {1., 0., 0., 0.};
  std::vector<double> p1 {0., 1., 0., 0.};
  std::vector<double> p2 {1., 1., 1., 0.};
  std::vector<double> p3 {1., 1., 0., 1.};
  std::vector<double> p4 {1., 1., -1., 0.};

  using Point_d = typename Unweighted_kernel::Point_d;
  std::vector<Point_d> unw_pts;
  unw_pts.emplace_back(p0.begin(), p0.end());
  unw_pts.emplace_back(p1.begin(), p1.end());
  unw_pts.emplace_back(p2.begin(), p2.end());
  unw_pts.emplace_back(p3.begin(), p3.end());
  unw_pts.emplace_back(p4.begin(), p4.end());

  Unweighted_kernel kernel;
  auto unw_sphere = kernel.get_sphere(unw_pts.cbegin(), unw_pts.cend());

  std::clog << "Center is " << unw_sphere.first << " - squared radius is " << unw_sphere.second << std::endl;

  using Weighted_kernel = Gudhi::alpha_complex::Alpha_kernel_d<CGAL_kernel, true>;

  using Weighted_point_d = typename Weighted_kernel::Weighted_point_d;
  using Bare_point_d = typename Weighted_kernel::Bare_point_d;
  std::vector<Weighted_point_d> w_pts;
  w_pts.emplace_back(Bare_point_d(p0.begin(), p0.end()), 0.);
  w_pts.emplace_back(Bare_point_d(p1.begin(), p1.end()), 0.);
  w_pts.emplace_back(Bare_point_d(p2.begin(), p2.end()), 0.);
  w_pts.emplace_back(Bare_point_d(p3.begin(), p3.end()), 0.);
  w_pts.emplace_back(Bare_point_d(p4.begin(), p4.end()), 0.);

  Weighted_kernel w_kernel;
  auto w_sphere = w_kernel.get_sphere(w_pts.cbegin(), w_pts.cend());

  std::clog << "Center is " << w_sphere.point() << " - squared radius is " << w_sphere.weight() << std::endl;

  CGAL::NT_converter<typename Weighted_kernel::FT, double> cast_to_double;
  // The results shall be the same with weights = 0.
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(cast_to_double(unw_sphere.second), cast_to_double(w_sphere.weight()));
  BOOST_CHECK(unw_sphere.first == w_sphere.point());

  auto unw_sq_rd = kernel.get_squared_radius(unw_pts.cbegin(), unw_pts.cend());
  std::clog << "Squared radius is " << unw_sq_rd << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(cast_to_double(unw_sphere.second), cast_to_double(unw_sq_rd));
  auto w_sq_rd = w_kernel.get_squared_radius(w_pts.cbegin(), w_pts.cend());
  std::clog << "Squared radius is " << w_sq_rd << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(cast_to_double(w_sphere.weight()), cast_to_double(w_sq_rd));
}
