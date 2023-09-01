/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "alpha_kernel_d_inexact_kernel_static"
#include <boost/test/unit_test.hpp>

#include <CGAL/Epick_d.h>

#include "Alpha_kernel_d_unit_test.h"

BOOST_AUTO_TEST_CASE(Alpha_kernel_d_inexact_kernel_static_dimension) {
  // test Alpha_kernel_d get_dimension with static dimension_tag
  test_alpha_kernel_d_dimension<CGAL::Epick_d< CGAL::Dimension_tag<4> >>();
}

BOOST_AUTO_TEST_CASE(Alpha_kernel_d_inexact_kernel_static_get_sphere) {
  // Use static dimension_tag for the user to be able to get sphere
  test_alpha_kernel_d_get_sphere<CGAL::Epick_d< CGAL::Dimension_tag<4> >>();
}
