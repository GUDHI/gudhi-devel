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
#define BOOST_TEST_MODULE "zero_weighted_alpha_complex_dynamic"
#include <boost/test/unit_test.hpp>

#include "Zero_weighted_alpha_complex_unit_test.h"

BOOST_AUTO_TEST_CASE(Zero_weighted_alpha_complex_dynamic) {
  do_test<CGAL::Epeck_d< CGAL::Dynamic_dimension_tag >>();
}
