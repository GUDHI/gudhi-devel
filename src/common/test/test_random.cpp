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
#include <numeric>  // for std::iota
#include <algorithm>  // for std::shuffle

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "random"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

using list_of_rnd_types = boost::mpl::list<double, float, int, unsigned int, long>;

BOOST_AUTO_TEST_CASE_TEMPLATE(random_same_seed, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_same_seed )" << std::endl;
  
  RndType min{0};
  RndType max{100};
  Gudhi::random::Random rng_1(1);
  RndType first  = rng_1.get(min, max);
  Gudhi::random::Random rng_2(1);
  RndType second = rng_2.get(min, max);

  std::clog << "First random number: " << first << " - Second random number: " << second << std::endl;
  BOOST_CHECK(first >= min);
  BOOST_CHECK(first <= max);
  BOOST_CHECK(first == second);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(random_different_seed, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_different_seed )" << std::endl;
  
  RndType min{0};
  RndType max{100};
  Gudhi::random::Random rng_1;
  RndType first  = rng_1.get(min, max);
  Gudhi::random::Random rng_2;
  RndType second = rng_2.get(min, max);

  std::clog << "First random number: " << first << " - Second random number: " << second << std::endl;
  BOOST_CHECK(first >= min);
  BOOST_CHECK(first <= max);
  BOOST_CHECK(second >= min);
  BOOST_CHECK(second <= max);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(default_random, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( default_random )" << std::endl;
  
  const RndType MIN_VAL{0};
  const RndType MAX_VAL{100};
  
  auto rng = Gudhi::random::get_default_random();
  RndType first  = rng.get(MIN_VAL, MAX_VAL);
  RndType second = rng.get(MIN_VAL, MAX_VAL);

  std::clog << "First random number: " << first << " - Second random number: " << second << std::endl;
  BOOST_CHECK(first >= MIN_VAL);
  BOOST_CHECK(first <= MAX_VAL);
  BOOST_CHECK(second >= MIN_VAL);
  BOOST_CHECK(second <= MAX_VAL);

  const int SEED{42};
  rng.set_seed(SEED);
  first  = rng.get(MIN_VAL, MAX_VAL);
  rng.set_seed(SEED);
  second = rng.get(MIN_VAL, MAX_VAL);

  std::clog << "First random number seeded : " << first << " - Second random number seeded : " << second << std::endl;
  BOOST_CHECK(first >= MIN_VAL);
  BOOST_CHECK(first <= MAX_VAL);
  BOOST_CHECK(first == second);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(operator_parenthesis, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( operator_parenthesis )" << std::endl;
  
  const RndType MIN_VAL{10};
  const RndType MAX_VAL{100};
  
  auto rng = Gudhi::random::get_default_random();
  RndType first  = rng.get(MIN_VAL, MAX_VAL);
  RndType second = rng.get(MIN_VAL, MAX_VAL);

  std::clog << ".operator(min, max) - First random number: " << first << " - Second random number: " << second << std::endl;
  BOOST_CHECK(first >= MIN_VAL);
  BOOST_CHECK(first <= MAX_VAL);
  BOOST_CHECK(second >= MIN_VAL);
  BOOST_CHECK(second <= MAX_VAL);

  first  = rng.get(MAX_VAL);
  second = rng.get(MAX_VAL);

  std::clog << ".operator() - First random number: " << first << " - Second random number: " << second << std::endl;
  BOOST_CHECK(first >= RndType{0});
  BOOST_CHECK(first <= MAX_VAL);
  BOOST_CHECK(second >= RndType{0});
  BOOST_CHECK(second <= MAX_VAL);

  first  = rng.get<RndType>();
  second = rng.get<RndType>();

  std::clog << ".operator() - First random number: " << first << " - Second random number: " << second << std::endl;
  BOOST_CHECK(first >= RndType{0});
  BOOST_CHECK(first <= RndType{1});
  BOOST_CHECK(second >= RndType{0});
  BOOST_CHECK(second <= RndType{1});
}

BOOST_AUTO_TEST_CASE(random_get_engine) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_get_engine )" << std::endl;
  
  std::vector<int> v(100);
  std::iota(v.begin(), v.end(), 0);
  auto rng = Gudhi::random::get_default_random();
  std::shuffle(v.begin(), v.end(), rng.get_engine());
  
  BOOST_CHECK(!std::is_sorted(v.cbegin(), v.cend()));
}
