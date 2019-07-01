/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carrière
 *
 *    Copyright (C) 2018  INRIA
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "kernel"

#include <boost/test/unit_test.hpp>
#include <cmath>  // float comparison
#include <limits>
#include <string>
#include <vector>
#include <algorithm>  // std::max
#include <gudhi/Persistence_heat_maps.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Sliced_Wasserstein.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>

using constant_scaling_function = Gudhi::Persistence_representations::constant_scaling_function;
using SW = Gudhi::Persistence_representations::Sliced_Wasserstein;
using PWG = Gudhi::Persistence_representations::Persistence_heat_maps<constant_scaling_function>;
using Persistence_diagram = std::vector<std::pair<double,double> >;

std::function<double(std::pair<double, double>, std::pair<double, double>)> Gaussian_function(double sigma){
  return [=](std::pair<double, double> p, std::pair<double, double> q){
    return (1/std::sqrt(2*Gudhi::Persistence_representations::pi)*sigma) * std::exp(  -( (p.first-q.first)*(p.first-q.first) + (p.second-q.second)*(p.second-q.second) )/(2*sigma)  );
  };
}

BOOST_AUTO_TEST_CASE(check_PWG) {
  Persistence_diagram v1, v2; v1.emplace_back(0,1); v2.emplace_back(0,2);
  PWG pwg1(v1, Gaussian_function(1.0));
  PWG pwg2(v2, Gaussian_function(1.0));
  BOOST_CHECK(std::abs(pwg1.compute_scalar_product(pwg2) - std::exp(-0.5)/(std::sqrt(2*Gudhi::Persistence_representations::pi))) <= 1e-3);
}

BOOST_AUTO_TEST_CASE(check_SW) {
  Persistence_diagram v1, v2; v1.emplace_back(0,1); v2.emplace_back(0,2);
  SW sw1(v1, 1.0, 100); SW swex1(v1, 1.0, -1);
  SW sw2(v2, 1.0, 100); SW swex2(v2, 1.0, -1);
  BOOST_CHECK(std::abs(sw1.compute_scalar_product(sw2) - swex1.compute_scalar_product(swex2)) <= 1e-1);
}
