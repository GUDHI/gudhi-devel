/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstddef>
#include <type_traits>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/Summand.h>
#include <gudhi/Multi_persistence/summand_helpers.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/Line.h>

using Gudhi::multi_persistence::Box;
using Gudhi::multi_persistence::Line;
using Gudhi::multi_persistence::Summand;
using Gudhi::multi_persistence::compute_summand_distance_to;
using Gudhi::multi_persistence::compute_summand_interleaving;
using Gudhi::multi_persistence::compute_summand_local_weight;
using Gudhi::multi_persistence::compute_summand_landscape_value;

template <typename T>
using B = typename Summand<T>::Births;
template <typename T>
using D = typename Summand<T>::Deaths;

using list_of_tested_variants = boost::mpl::list<double, float, int, unsigned int>;

BOOST_AUTO_TEST_CASE_TEMPLATE(summand_constructors, T, list_of_tested_variants) {
  const int numParam = 3;

  Summand<T> empty(numParam);
  BOOST_CHECK_EQUAL(empty.get_dimension(), Summand<T>::template get_null_value<int>());
  BOOST_CHECK_EQUAL(empty.get_number_of_parameters(), numParam);
  BOOST_CHECK_EQUAL(empty.get_number_of_birth_corners(), 1);
  BOOST_CHECK_EQUAL(empty.get_number_of_death_corners(), 1);
  BOOST_CHECK(empty.get_upset() == B<T>::inf(numParam));
  BOOST_CHECK(empty.get_downset() == D<T>::minus_inf(numParam));

  std::vector<T> b = {0, 1, 2, 3, 4, 5};
  std::vector<T> d = {6, 7, 8};

  Summand<T> sum1(B<T>(b.begin(), b.end(), numParam), D<T>(d.begin(), d.end(), numParam), 2);
  BOOST_CHECK_EQUAL(sum1.get_dimension(), 2);
  BOOST_CHECK_EQUAL(sum1.get_number_of_parameters(), numParam);
  BOOST_CHECK_EQUAL(sum1.get_number_of_birth_corners(), 2);
  BOOST_CHECK_EQUAL(sum1.get_number_of_death_corners(), 1);
  const auto& upset1 = sum1.get_upset();
  const auto& downset1 = sum1.get_downset();
  BOOST_CHECK_EQUAL(upset1(0, 0), 0);
  BOOST_CHECK_EQUAL(upset1(0, 1), 1);
  BOOST_CHECK_EQUAL(upset1(0, 2), 2);
  BOOST_CHECK_EQUAL(upset1(1, 0), 3);
  BOOST_CHECK_EQUAL(upset1(1, 1), 4);
  BOOST_CHECK_EQUAL(upset1(1, 2), 5);
  BOOST_CHECK_EQUAL(downset1(0, 0), 6);
  BOOST_CHECK_EQUAL(downset1(0, 1), 7);
  BOOST_CHECK_EQUAL(downset1(0, 2), 8);

  const auto& flatUpset = sum1.compute_flat_upset();
  const auto& flatDownset = sum1.compute_flat_downset();
  BOOST_CHECK(flatUpset == b);
  BOOST_CHECK(flatDownset == d);

  Summand<T> sum2(b, d, numParam, 2);
  BOOST_CHECK_EQUAL(sum2.get_dimension(), 2);
  BOOST_CHECK_EQUAL(sum2.get_number_of_parameters(), numParam);
  BOOST_CHECK_EQUAL(sum2.get_number_of_birth_corners(), 2);
  BOOST_CHECK_EQUAL(sum2.get_number_of_death_corners(), 1);
  const auto& upset2 = sum2.get_upset();
  const auto& downset2 = sum2.get_downset();
  BOOST_CHECK_EQUAL(upset2(0, 0), 0);
  BOOST_CHECK_EQUAL(upset2(0, 1), 1);
  BOOST_CHECK_EQUAL(upset2(0, 2), 2);
  BOOST_CHECK_EQUAL(upset2(1, 0), 3);
  BOOST_CHECK_EQUAL(upset2(1, 1), 4);
  BOOST_CHECK_EQUAL(upset2(1, 2), 5);
  BOOST_CHECK_EQUAL(downset2(0, 0), 6);
  BOOST_CHECK_EQUAL(downset2(0, 1), 7);
  BOOST_CHECK_EQUAL(downset2(0, 2), 8);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(summand_access, T, list_of_tested_variants) {
  using P = typename Box<T>::Point_t;

  std::vector<T> b = {1, 3, 3, 1};
  std::vector<T> d = {6, 6};

  Summand<T> sum(b, d, 2, 2);
  Box<T> box = sum.compute_bounds();
  BOOST_CHECK((box.get_lower_corner() == P{1, 1}));
  BOOST_CHECK((box.get_upper_corner() == P{6, 6}));

  BOOST_CHECK(sum.contains(B<T>({2, 5})));
  BOOST_CHECK(sum.contains(B<T>({5, 2})));
  BOOST_CHECK(sum.contains(B<T>({4, 4})));
  BOOST_CHECK(!sum.contains(B<T>({0, 5})));
  BOOST_CHECK(!sum.contains(B<T>({1, 1})));
  BOOST_CHECK(!sum.contains(B<T>({8, 5})));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(summand_corner_modifiers, T, list_of_tested_variants) {
  std::vector<T> b = {2, 3, 4, 1};
  std::vector<T> d = {7, 6};

  Summand<T> sum(b, d, 2, 2);

  const auto& upset = sum.get_upset();
  const auto& downset = sum.get_downset();
  BOOST_CHECK_EQUAL(upset.num_generators(), 2);
  BOOST_CHECK_EQUAL(upset(0, 0), 2);
  BOOST_CHECK_EQUAL(upset(0, 1), 3);
  BOOST_CHECK_EQUAL(upset(1, 0), 4);
  BOOST_CHECK_EQUAL(upset(1, 1), 1);
  BOOST_CHECK_EQUAL(downset.num_generators(), 1);
  BOOST_CHECK_EQUAL(downset(0, 0), 7);
  BOOST_CHECK_EQUAL(downset(0, 1), 6);

  sum.add_bar({6, 0}, {9, 4});
  BOOST_CHECK_EQUAL(upset.num_generators(), 3);
  BOOST_CHECK_EQUAL(upset(0, 0), 2);
  BOOST_CHECK_EQUAL(upset(0, 1), 3);
  BOOST_CHECK_EQUAL(upset(1, 0), 4);
  BOOST_CHECK_EQUAL(upset(1, 1), 1);
  BOOST_CHECK_EQUAL(upset(2, 0), 6);
  BOOST_CHECK_EQUAL(upset(2, 1), 0);
  BOOST_CHECK_EQUAL(downset.num_generators(), 2);
  BOOST_CHECK_EQUAL(downset(0, 0), 7);
  BOOST_CHECK_EQUAL(downset(0, 1), 6);
  BOOST_CHECK_EQUAL(downset(1, 0), 9);
  BOOST_CHECK_EQUAL(downset(1, 1), 4);

  Line<T> l({2, 6}, {2, 1});
  Box<T> box({1, 4}, {5, 7});

  auto bar = sum.get_bar(Line<T>({2, 7}, {2, 1}));
  BOOST_CHECK_EQUAL(bar[0], std::numeric_limits<double>::infinity());
  BOOST_CHECK_EQUAL(bar[1], std::numeric_limits<double>::infinity());
  bar = sum.get_bar(l);
  BOOST_CHECK_EQUAL(bar[0], 0);
  BOOST_CHECK_EQUAL(bar[1], 0);

  sum.add_bar(l, -1, 10, box, true);
  BOOST_CHECK_EQUAL(upset.num_generators(), 4);
  BOOST_CHECK_EQUAL(upset(0, 0), 1);
  BOOST_CHECK_EQUAL(upset(0, 1), 5);
  BOOST_CHECK_EQUAL(upset(1, 0), 2);
  BOOST_CHECK_EQUAL(upset(1, 1), 3);
  BOOST_CHECK_EQUAL(upset(2, 0), 4);
  BOOST_CHECK_EQUAL(upset(2, 1), 1);
  BOOST_CHECK_EQUAL(upset(3, 0), 6);
  BOOST_CHECK_EQUAL(upset(3, 1), 0);
  BOOST_CHECK_EQUAL(downset.num_generators(), 3);
  BOOST_CHECK_EQUAL(downset(0, 0), 5);
  BOOST_CHECK_EQUAL(downset(0, 1), 7);
  BOOST_CHECK_EQUAL(downset(1, 0), 7);
  BOOST_CHECK_EQUAL(downset(1, 1), 6);
  BOOST_CHECK_EQUAL(downset(2, 0), 9);
  BOOST_CHECK_EQUAL(downset(2, 1), 4);

  bar = sum.get_bar(l);
  BOOST_CHECK_EQUAL(bar[0], -0.5);
  BOOST_CHECK_EQUAL(bar[1], 1.);

  sum.identify_births(1);
  BOOST_CHECK_EQUAL(upset.num_generators(), 4);
  BOOST_CHECK_EQUAL(upset(0, 0), 1);
  BOOST_CHECK_EQUAL(upset(0, 1), 5);
  BOOST_CHECK_EQUAL(upset(1, 0), 2);
  BOOST_CHECK_EQUAL(upset(1, 1), 3);
  BOOST_CHECK_EQUAL(upset(2, 0), 4);
  BOOST_CHECK_EQUAL(upset(2, 1), 1);
  BOOST_CHECK_EQUAL(upset(3, 0), 6);
  BOOST_CHECK_EQUAL(upset(3, 1), 0);
  sum.identify_deaths(1);
  BOOST_CHECK_EQUAL(downset.num_generators(), 3);
  BOOST_CHECK_EQUAL(downset(0, 0), 5);
  BOOST_CHECK_EQUAL(downset(0, 1), 7);
  BOOST_CHECK_EQUAL(downset(1, 0), 7);
  BOOST_CHECK_EQUAL(downset(1, 1), 6);
  BOOST_CHECK_EQUAL(downset(2, 0), 9);
  BOOST_CHECK_EQUAL(downset(2, 1), 4);

  sum.identify_births(3);
  BOOST_CHECK_EQUAL(upset.num_generators(), 2);
  BOOST_CHECK_EQUAL(upset(0, 0), 1);
  BOOST_CHECK_EQUAL(upset(0, 1), 3);
  BOOST_CHECK_EQUAL(upset(1, 0), 4);
  BOOST_CHECK_EQUAL(upset(1, 1), 0);
  sum.identify_deaths(2);
  BOOST_CHECK_EQUAL(downset.num_generators(), 3);
  BOOST_CHECK_EQUAL(downset(0, 0), 5);
  BOOST_CHECK_EQUAL(downset(0, 1), 7);
  BOOST_CHECK_EQUAL(downset(1, 0), 7);
  BOOST_CHECK_EQUAL(downset(1, 1), 6);
  BOOST_CHECK_EQUAL(downset(2, 0), 9);
  BOOST_CHECK_EQUAL(downset(2, 1), 4);

  sum.add_bar(l, -1, 10, box, false);
  BOOST_CHECK_EQUAL(upset.num_generators(), 3);
  BOOST_CHECK_EQUAL(upset(0, 0), Summand<T>::T_m_inf);
  BOOST_CHECK_EQUAL(upset(0, 1), 5);
  BOOST_CHECK_EQUAL(upset(1, 0), 1);
  BOOST_CHECK_EQUAL(upset(1, 1), 3);
  BOOST_CHECK_EQUAL(upset(2, 0), 4);
  BOOST_CHECK_EQUAL(upset(2, 1), 0);
  BOOST_CHECK_EQUAL(downset.num_generators(), 1);
  BOOST_CHECK(downset == D<T>::inf(2));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(summand_transformers, T, list_of_tested_variants) {
  std::vector<T> b = {2, 3, 4, 1};
  std::vector<T> d = {7, 6};

  Summand<T> sum(b, d, 2, 2);
  std::vector<std::vector<T>> grid = {{2, 5, 6, 9, 11, 12, 16}, {3, 4, 6, 7, 10, 15, 20, 22, 26}};
  sum.evaluate_in_grid(grid);

  const auto& upset = sum.get_upset();
  const auto& downset = sum.get_downset();
  BOOST_CHECK_EQUAL(upset.num_generators(), 2);
  BOOST_CHECK_EQUAL(upset(0, 0), 6);
  BOOST_CHECK_EQUAL(upset(0, 1), 7);
  BOOST_CHECK_EQUAL(upset(1, 0), 11);
  BOOST_CHECK_EQUAL(upset(1, 1), 4);
  BOOST_CHECK_EQUAL(downset.num_generators(), 1);
  BOOST_CHECK_EQUAL(downset(0, 0), Summand<T>::T_inf);
  BOOST_CHECK_EQUAL(downset(0, 1), 20);

  sum.rescale(std::vector<double>{2, 3});
  BOOST_CHECK_EQUAL(upset.num_generators(), 2);
  BOOST_CHECK_EQUAL(upset(0, 0), 12);
  BOOST_CHECK_EQUAL(upset(0, 1), 21);
  BOOST_CHECK_EQUAL(upset(1, 0), 22);
  BOOST_CHECK_EQUAL(upset(1, 1), 12);
  BOOST_CHECK_EQUAL(downset.num_generators(), 1);
  BOOST_CHECK_EQUAL(downset(0, 0), Summand<T>::T_inf);
  BOOST_CHECK_EQUAL(downset(0, 1), 60);

  sum.translate(std::vector<double>{-10, 3});
  BOOST_CHECK_EQUAL(upset.num_generators(), 2);
  BOOST_CHECK_EQUAL(upset(0, 0), 2);
  BOOST_CHECK_EQUAL(upset(0, 1), 24);
  BOOST_CHECK_EQUAL(upset(1, 0), 12);
  BOOST_CHECK_EQUAL(upset(1, 1), 15);
  BOOST_CHECK_EQUAL(downset.num_generators(), 1);
  BOOST_CHECK_EQUAL(downset(0, 0), Summand<T>::T_inf);
  BOOST_CHECK_EQUAL(downset(0, 1), 63);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(summand_helpers, T, list_of_tested_variants) {
  std::vector<T> b = {2, 3, 4, 1};
  std::vector<T> d = {3, 8, 7, 6};

  Summand<T> sum(b, d, 2, 2);
  auto distance = compute_summand_distance_to(sum, std::vector<T>{1, 4}, false);
  BOOST_CHECK_EQUAL(distance, 1);
  distance = compute_summand_distance_to(sum, std::vector<T>{5, 5}, false);
  BOOST_CHECK_EQUAL(distance, 0);
  distance = compute_summand_distance_to(sum, std::vector<T>{10, 3}, false);
  BOOST_CHECK_EQUAL(distance, 3);
  distance = compute_summand_distance_to(sum, std::vector<T>{1, 4}, true);
  BOOST_CHECK_EQUAL(distance, 1);
  distance = compute_summand_distance_to(sum, std::vector<T>{5, 5}, true);
  BOOST_CHECK_EQUAL(distance, -1);
  distance = compute_summand_distance_to(sum, std::vector<T>{10, 3}, true);
  BOOST_CHECK_EQUAL(distance, 3);

  auto inter = compute_summand_interleaving(sum, Box<T>());
  BOOST_CHECK_EQUAL(inter, 3);
  inter = compute_summand_interleaving(sum, Box<T>({1, 1}, {7, 8}));
  BOOST_CHECK_EQUAL(inter, 3);
  inter = compute_summand_interleaving(sum, Box<T>({3, 2}, {6, 5}));
  BOOST_CHECK_EQUAL(inter, 2);

  auto weight = compute_summand_local_weight(sum, std::vector<T>{4, 3}, 4.);
  BOOST_CHECK_EQUAL(weight, 0.375);
  weight = compute_summand_local_weight(sum, std::vector<T>{5, 2}, 2.);
  BOOST_CHECK_EQUAL(weight, 0.75);
  if constexpr (std::is_integral_v<T>) {
    weight = compute_summand_local_weight(sum, std::vector<T>{5, 2}, 0.75);
    BOOST_CHECK_EQUAL(weight, 1. / 0.75);
  } else {
    weight = compute_summand_local_weight(sum, std::vector<T>{5, 2}, 0.75);
    BOOST_CHECK_EQUAL(weight, 1.);
  }
  weight = compute_summand_local_weight(sum, std::vector<T>{4, 3}, -4.);
  BOOST_CHECK_EQUAL(weight, 15. / 64.);
  weight = compute_summand_local_weight(sum, std::vector<T>{5, 2}, -2.);
  BOOST_CHECK_EQUAL(weight, 9. / 16.);
  if constexpr (std::is_integral_v<T>) {
    weight = compute_summand_local_weight(sum, std::vector<T>{5, 2}, -0.75);
    BOOST_CHECK_EQUAL(weight, 4. / 2.25);
  } else {
    weight = compute_summand_local_weight(sum, std::vector<T>{5, 2}, -0.75);
    BOOST_CHECK_EQUAL(weight, 1.);
  }

  auto landscape = compute_summand_landscape_value(sum, std::vector<T>{1, 4});
  BOOST_CHECK_EQUAL(landscape, 0);
  landscape = compute_summand_landscape_value(sum, std::vector<T>{5, 5});
  BOOST_CHECK_EQUAL(landscape, 1);
  landscape = compute_summand_landscape_value(sum, std::vector<T>{10, 3});
  BOOST_CHECK_EQUAL(landscape, 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(summand_serialization, T, list_of_tested_variants) {
  std::vector<T> b = {2, 3, 4, 1};
  std::vector<T> d = {3, 8, 7, 6, 6, 99};

  Summand<T> sum(b, d, 2, 2);
  char* buffer = new char[256];
  char* ptr = buffer;

  std::size_t serSize = get_serialization_size_of(sum);
  ptr = serialize_value_to_char_buffer(sum, ptr);
  BOOST_CHECK_EQUAL(serSize, ptr - buffer);

  Summand<T> copy;
  const char* c_ptr = buffer;
  c_ptr = deserialize_value_from_char_buffer(copy, c_ptr);
  BOOST_CHECK_EQUAL(serSize, c_ptr - buffer);
  BOOST_CHECK(sum == copy);

  delete [] buffer;
}
