/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA (France)
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
#define BOOST_TEST_MODULE "Persistence_representations"
#include <boost/test/unit_test.hpp>
#include <gudhi/reader_utils.h>
#include "gudhi/Persistence_intervals_with_distances.h"
#include "gudhi/common_persistence_representations.h"

#include <iostream>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

BOOST_AUTO_TEST_CASE(check_bottleneck_distances_computation) {
  Persistence_intervals_with_distances p("data/file_with_diagram");
  Persistence_intervals_with_distances q("data/file_with_diagram_1");

  double dist = p.distance(q);

  // std::cout << "dist : " << dist << std::endl;

  BOOST_CHECK(almost_equal(dist, 0.389043));
}

BOOST_AUTO_TEST_CASE(check_default_parameters_in_distance) {
  Persistence_intervals_with_distances p("data/file_with_diagram");
  Persistence_intervals_with_distances q("data/file_with_diagram_1");

  double default_parameter_distance = p.distance(q);
  double max_parameter_distance = p.distance(q, std::numeric_limits<double>::max());
  double inf_parameter_distance = p.distance(q, std::numeric_limits<double>::infinity());

  // std::cout << "default_parameter_distance : " << default_parameter_distance << std::endl;
  // std::cout << "max_parameter_distance : " << max_parameter_distance << std::endl;
  // std::cout << "inf_parameter_distance : " << inf_parameter_distance << std::endl;

  BOOST_CHECK(default_parameter_distance == max_parameter_distance);
  BOOST_CHECK(inf_parameter_distance == max_parameter_distance);
  BOOST_CHECK(inf_parameter_distance == max_parameter_distance);
}
