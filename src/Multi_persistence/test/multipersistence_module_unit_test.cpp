/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

// #include <array>
// #include <initializer_list>
// #include <limits>
// #include <utility>
// #include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/Module.h>

using Gudhi::multi_persistence::Module;

using list_of_tested_variants = boost::mpl::list<double, float, int, unsigned int>;

BOOST_AUTO_TEST_CASE_TEMPLATE(module_constructors, T, list_of_tested_variants)
{
  Module<T> empty;
  BOOST_CHECK_EQUAL(empty.size(), 0);
  BOOST_CHECK(empty.begin() == empty.end());
}

