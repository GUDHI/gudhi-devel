/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "delaunay_complex_exact_kernel_static"
#include <boost/test/unit_test.hpp>

#include <CGAL/Epeck_d.h>

#include "Delaunay_complex_unit_test.h"

BOOST_AUTO_TEST_CASE(Delaunay_complex_exact_kernel_static_simplices_comparison) {
  // Use static dimension_tag for the user not to be able to set dimension
  compare_delaunay_complex_simplices<CGAL::Epeck_d< CGAL::Dimension_tag<5> >>();
}
