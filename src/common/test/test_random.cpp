/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Random.h>

#include <iostream>
#include <vector>
#include <numeric>  // for std::iota
#include <algorithm>  // for std::shuffle
#include <random>  // for std::random_device

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "random"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

using list_of_rnd_types = boost::mpl::list<double, float, int, unsigned int, long>;

// Note: If you need to add some tests, consider to separate them, with first, all tests that are requiring some
// randomness, and then, the tests that sets the seed (cf. TESTS_WITH_SEED).

// TESTS_WITHOUT_SEED
BOOST_AUTO_TEST_CASE_TEMPLATE(random_get_limits, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_get_limits )" << "\n";
  
  const RndType MIN_VAL{10};
  const RndType MAX_VAL{100};
  
  RndType random_value = Gudhi::random::get_uniform(MIN_VAL, MAX_VAL);

  std::clog << "get_default_random().get(min, max) - random value: " << random_value << "\n";
  BOOST_CHECK(random_value >= MIN_VAL);
  BOOST_CHECK(random_value <= MAX_VAL);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(random_get_range_limits, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_get_limits )" << "\n";
  
  const RndType MIN_VAL{10};
  const RndType MAX_VAL{100};
  int NB{50};
  
  std::vector<RndType> first_range   = Gudhi::random::get_uniform_range(NB, MIN_VAL, MAX_VAL);
  std::vector<RndType> second_range  = Gudhi::random::get_uniform_range(NB, MIN_VAL, MAX_VAL);

  std::clog << "get_default_random().get_range(nb, min, max)\n";
  for (auto val: first_range) {
    std::clog << val << ", ";
    BOOST_CHECK(val >= MIN_VAL);
    BOOST_CHECK(val <= MAX_VAL);
  }
  std::clog << "\n";
  for (auto val: second_range) {
    std::clog << val << ", ";
    BOOST_CHECK(val >= MIN_VAL);
    BOOST_CHECK(val <= MAX_VAL);
  }
  std::clog << "\n";
  // NB times the same values is not normal
  BOOST_CHECK(first_range != second_range);
}

BOOST_AUTO_TEST_CASE(random_generator_external_use) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_generator_external_use )" << "\n";
  
  std::vector<int> v(100);
  std::iota(v.begin(), v.end(), 0);
  std::shuffle(v.begin(), v.end(), Gudhi::random::get_default_random());
  
  BOOST_CHECK(!std::is_sorted(v.cbegin(), v.cend()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(custom_random_generator, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( custom_random_generator )" << "\n";
  const RndType MIN{5};
  const RndType MAX{10};
  const int NB{50};
  auto rng = std::default_random_engine{std::random_device{}()};
  std::clog << "Custom random engine - get() returns " << Gudhi::random::get_uniform(MIN, MAX, rng) << "\n";
  
  auto first_range  = Gudhi::random::get_uniform_range(NB, MIN, MAX, rng);
  auto second_range = Gudhi::random::get_uniform_range(NB, MIN, MAX, rng);
  
  std::clog << "Custom random engine - get_range()\nFirst range:\n";
  for (RndType value: first_range) {
    std::clog << value << ", ";
  }
  std::clog << "\nSecond range:\n";
  for (RndType value: second_range) {
    std::clog << value << ", ";
  }
  std::clog << "\n";
  // No seed means at least one value on the NB ones should be different
  BOOST_CHECK(first_range != second_range);
  
  rng.seed(42);
  auto first  = Gudhi::random::get_uniform(MIN, MAX, rng);
  rng.seed(42);
  auto second = Gudhi::random::get_uniform(MIN, MAX, rng);
  std::clog << "Custom random engine with seed - first = " << first << " - second = " << second << "\n";
  BOOST_CHECK(first == second);

  rng.seed(42);
  first_range  = Gudhi::random::get_uniform_range(NB, MIN, MAX, rng);
  rng.seed(42);
  second_range = Gudhi::random::get_uniform_range(NB, MIN, MAX, rng);
  
  std::clog << "Custom random engine - get_range() with a seed\nFirst range:\n";
  for (RndType value: first_range) {
    std::clog << value << ", ";
  }
  std::clog << "\nSecond range:\n";
  for (RndType value: second_range) {
    std::clog << value << ", ";
  }
  std::clog << "\n";
  BOOST_CHECK(first_range == second_range);
}

// TESTS_WITH_SEED

BOOST_AUTO_TEST_CASE_TEMPLATE(random_get_with_seed, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_get_with_seed )" << "\n";
  
  const RndType MIN_VAL{10};
  const RndType MAX_VAL{100};
  Gudhi::random::set_seed(42);
  RndType first  = Gudhi::random::get_uniform(MIN_VAL, MAX_VAL);
  Gudhi::random::set_seed(42);
  RndType second = Gudhi::random::get_uniform(MIN_VAL, MAX_VAL);

  std::clog << "Gudhi::random::set_seed + Gudhi::random::get_uniform(min, max) - First random number: " << first
            << " - Second random number: " << second << "\n";
  BOOST_CHECK(first >= MIN_VAL);
  BOOST_CHECK(first <= MAX_VAL);
  BOOST_CHECK(first == second);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(random_get_range_with_seed, RndType, list_of_rnd_types) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_get_range_with_seed )" << "\n";
  
  const RndType MIN_VAL{10};
  const RndType MAX_VAL{100};
  int NB{50};
  Gudhi::random::set_seed(42);
  std::vector<RndType> first_range  = Gudhi::random::get_uniform_range(NB, MIN_VAL, MAX_VAL);
  Gudhi::random::set_seed(42);
  std::vector<RndType> second_range = Gudhi::random::get_uniform_range(NB, MIN_VAL, MAX_VAL);

  std::clog << "Gudhi::random::set_seed + Gudhi::random::get_uniform_range(nb, min, max)\n";
  for (auto val: first_range) {
    std::clog << val << ", ";
    BOOST_CHECK(val >= MIN_VAL);
    BOOST_CHECK(val <= MAX_VAL);
  }
  std::clog << "\n";
  for (auto val: second_range) {
    std::clog << val << ", ";
    BOOST_CHECK(val >= MIN_VAL);
    BOOST_CHECK(val <= MAX_VAL);
  }
  std::clog << "\n";
  BOOST_CHECK(first_range == second_range);
}

BOOST_AUTO_TEST_CASE(random_generator_external_use_with_seed) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_generator_external_use_with_seed )" << "\n";
  
  Gudhi::random::set_seed(42);
  std::vector<int> vec1(100);
  std::iota(vec1.begin(), vec1.end(), 0);
  std::shuffle(vec1.begin(), vec1.end(), Gudhi::random::get_default_random());

  Gudhi::random::set_seed(42);
  std::vector<int> vec2(100);
  std::iota(vec2.begin(), vec2.end(), 0);
  std::shuffle(vec2.begin(), vec2.end(), Gudhi::random::get_default_random());
  
  BOOST_CHECK(vec1 == vec2);
}
