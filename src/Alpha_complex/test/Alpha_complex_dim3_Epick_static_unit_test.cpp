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
#define BOOST_TEST_MODULE "alpha_complex_dim3_inexact_kernel_static"
#include <boost/test/unit_test.hpp>

#include <CGAL/Epick_d.h>

#include "Alpha_complex_dim3_unit_test.h"

// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dimension_tag<3> > Inexact_kernel_s;

BOOST_AUTO_TEST_CASE(Alpha_complex_from_OFF_file_inexact_kernel_static_dimension) {
  test_alpha_complex_from_OFF_file<Inexact_kernel_s>();
}

BOOST_AUTO_TEST_CASE(Alpha_complex_from_empty_points_inexact_kernel_static_dimension) {
  test_alpha_complex_from_empty_points<Inexact_kernel_s>();
}
