/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>  // std::pair, std::make_pair
#include <cmath>  // float comparison
#include <limits>
#include <functional>  // greater
#include <tuple>  // std::tie

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree"
#include <boost/test/unit_test.hpp>

#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

BOOST_AUTO_TEST_CASE(euler) {
  typedef Simplex_tree<> ST;
  ST st;
  BOOST_CHECK(st.euler_characteristic() == 0);
  BOOST_CHECK(st.magnitude(-.5) == 0.);
  st.insert_simplex_and_subfaces({0,1,2},.5);
  BOOST_CHECK(st.euler_characteristic() == 1);
  BOOST_CHECK(std::abs(st.magnitude(0) - 1.) < .000001);
  BOOST_CHECK(std::abs(st.magnitude(2.) - std::exp(-1.)) < .000001);
  st.insert_simplex({0},0);
  BOOST_CHECK(st.euler_characteristic() == 1);
  BOOST_CHECK(std::abs(st.magnitude(2.) - 1.) < .000001);
}
