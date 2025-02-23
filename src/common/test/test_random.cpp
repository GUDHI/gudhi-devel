/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Random.h>

#include <iostream>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "random"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

using list_of_rnd_types = boost::mpl::list<double, float, int, unsigned int, long>;

BOOST_AUTO_TEST_CASE_TEMPLATE(random_same_seed, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_same_seed )" << std::endl;
  
  RndType min{0};
  RndType max{100};
  RndType first  = Gudhi::Random(1).get<RndType>(min, max);
  RndType second = Gudhi::Random(1).get<RndType>(min, max);

  std::clog << "First random number: " << first << " - Second random number: " << second << std::endl;
  BOOST_CHECK(first >= min);
  BOOST_CHECK(first <= max);
  BOOST_CHECK(first == second);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(random_different_seed, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_different_seed )" << std::endl;
  
  RndType min{0};
  RndType max{100};
  RndType first  = Gudhi::Random().get<RndType>(min, max);
  RndType second = Gudhi::Random().get<RndType>(min, max);

  std::clog << "First random number: " << first << " - Second random number: " << second << std::endl;
  BOOST_CHECK(first >= min);
  BOOST_CHECK(first <= max);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(random_range, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_range )" << std::endl;
  
  RndType min{0};
  RndType max{100};
  const std::size_t N{20};
  std::vector<RndType> range = Gudhi::Random().get_range<RndType>(N, min, max);
  
  BOOST_CHECK(range.size() == N);

  for (const auto& value : range) {
    std::clog << "Random number: " << value << std::endl;
    BOOST_CHECK(value >= min);
    BOOST_CHECK(value <= max);
  }
}
