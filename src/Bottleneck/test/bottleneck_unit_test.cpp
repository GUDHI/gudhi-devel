#define BOOST_TEST_MODULE bottleneck test

#include <boost/test/included/unit_test.hpp>

#include "gudhi/Graph_matching.h"
#include <iostream>

using namespace Gudhi::bottleneck;

BOOST_AUTO_TEST_CASE(random_diagrams) {
  int n = 100;
  // Random construction
  std::vector< std::pair<double, double> > v1, v2;
  for (int i = 0; i < n; i++) {
    int a = rand() % n;
    v1.emplace_back(a, a + rand() % (n - a));
    int b = rand() % n;
    v2.emplace_back(b, b + rand() % (n - b));
  }
  // v1 and v2 are persistence diagrams containing each 100 randoms points.
  double b = bottleneck_distance(v1, v2, 0);
  //
  std::cout << b << std::endl;
  const double EXPECTED_DISTANCE = 98.5;
  BOOST_CHECK(b == EXPECTED_DISTANCE);
}
