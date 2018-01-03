/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carrière
 *
 *    Copyright (C) 2017  INRIA
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "kernel"

#include <boost/test/unit_test.hpp>
#include <cmath>  // float comparison
#include <limits>
#include <string>
#include <vector>
#include <algorithm>  // std::max
#include <gudhi/kernel.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>

BOOST_AUTO_TEST_CASE(check_PSS) {
  std::vector< std::pair<double, double> > v1, v2;
  v1.emplace_back(std::pair<double,double>(0,1));
  v2.emplace_back(std::pair<double,double>(0,2));
  BOOST_CHECK(std::abs(Gudhi::kernel::pssk(v1,v2,1) - Gudhi::kernel::approx_pssk(v1,v2,1)) <= 1e-1);
}

BOOST_AUTO_TEST_CASE(check_PWG) {
  std::vector< std::pair<double, double> > v1, v2;
  v1.emplace_back(std::pair<double,double>(0,1));
  v2.emplace_back(std::pair<double,double>(0,2));
  BOOST_CHECK(std::abs(Gudhi::kernel::lpwgk(v1,v2,1) - Gudhi::kernel::approx_lpwgk(v1,v2,1)) <= 1e-1);
  BOOST_CHECK(std::abs(Gudhi::kernel::gpwgk(v1,v2,1,1) - Gudhi::kernel::approx_gpwgk(v1,v2,1,1)) <= 1e-1);
}

BOOST_AUTO_TEST_CASE(check_SW) {
  std::vector< std::pair<double, double> > v1, v2;
  v2.emplace_back(std::pair<double,double>(0,2));
  BOOST_CHECK(std::abs(Gudhi::kernel::sw(v1,v2) - Gudhi::kernel::approx_sw(v1,v2)) <= 1e-3);
  BOOST_CHECK(std::abs(Gudhi::kernel::sw(v1,v2) - 2*std::sqrt(2)/3.1415) <= 1e-3);
}
