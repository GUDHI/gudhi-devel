/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Persistence_intervals_with_distances_test"
#include <boost/test/unit_test.hpp>
#include <gudhi/reader_utils.h>
#include <gudhi/Persistence_intervals_with_distances.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Unitary_tests_utils.h>

#include <iostream>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

BOOST_AUTO_TEST_CASE(check_bottleneck_distances_computation) {
  Persistence_intervals_with_distances p("data/file_with_diagram");
  Persistence_intervals_with_distances q("data/file_with_diagram_1");

  double dist = p.distance(q);

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(dist, 0.389043, Gudhi::Persistence_representations::epsi);
}

BOOST_AUTO_TEST_CASE(check_default_parameters_in_distance) {
  Persistence_intervals_with_distances p("data/file_with_diagram");
  Persistence_intervals_with_distances q("data/file_with_diagram_1");

  double default_parameter_distance = p.distance(q);
  double max_parameter_distance = p.distance(q, std::numeric_limits<double>::max());
  double inf_parameter_distance = p.distance(q, std::numeric_limits<double>::infinity());

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(default_parameter_distance, max_parameter_distance);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(inf_parameter_distance, max_parameter_distance);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(inf_parameter_distance, max_parameter_distance);
}
