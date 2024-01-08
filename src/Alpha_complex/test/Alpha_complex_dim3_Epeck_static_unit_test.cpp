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
#define BOOST_TEST_MODULE "alpha_complex_dim3_exact_kernel_static"
#include <boost/test/unit_test.hpp>

#include <CGAL/Epeck_d.h>

#include "Alpha_complex_dim3_unit_test.h"

// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epeck_d< CGAL::Dimension_tag<3> > Exact_kernel_s;

BOOST_AUTO_TEST_CASE(Alpha_complex_from_OFF_file_exact_kernel_static_dimension) {
  test_alpha_complex_from_OFF_file<Exact_kernel_s>();
}

BOOST_AUTO_TEST_CASE(Alpha_complex_from_empty_points_exact_kernel_static_dimension) {
  test_alpha_complex_from_empty_points<Exact_kernel_s>();
}
