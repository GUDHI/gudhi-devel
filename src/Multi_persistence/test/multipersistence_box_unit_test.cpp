/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024-25 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <limits>
#include <utility>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/Box.h>

using Gudhi::multi_persistence::Box;

typedef boost::mpl::list<double, float, int> list_of_tested_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(box_constructors, T, list_of_tested_variants)
{
  Box<T> b;
  BOOST_CHECK(b.is_trivial());
  auto bottom = b.get_lower_corner();
  auto top = b.get_upper_corner();
  BOOST_CHECK_EQUAL(bottom.size(), 0);
  BOOST_CHECK_EQUAL(top.size(), 0);
  BOOST_CHECK(b.get_bounding_corners().first == bottom);
  BOOST_CHECK(b.get_bounding_corners().second == top);

  Box<T> b2({0, 1, 2}, {4, 5, 6});
  BOOST_CHECK(!b2.is_trivial());
  bottom = b2.get_lower_corner();
  top = b2.get_upper_corner();
  BOOST_CHECK_EQUAL(bottom.size(), 3);
  BOOST_CHECK_EQUAL(top.size(), 3);
  BOOST_CHECK_EQUAL(bottom[0], 0);
  BOOST_CHECK_EQUAL(bottom[1], 1);
  BOOST_CHECK_EQUAL(bottom[2], 2);
  BOOST_CHECK_EQUAL(top[0], 4);
  BOOST_CHECK_EQUAL(top[1], 5);
  BOOST_CHECK_EQUAL(top[2], 6);
  BOOST_CHECK(b2.get_bounding_corners().first == bottom);
  BOOST_CHECK(b2.get_bounding_corners().second == top);

  Box<T> b3(std::make_pair<typename Box<T>::Point_t, typename Box<T>::Point_t>({0, 1, 2}, {4, 5, 6}));
  BOOST_CHECK(!b3.is_trivial());
  bottom = b3.get_lower_corner();
  top = b3.get_upper_corner();
  BOOST_CHECK_EQUAL(bottom.size(), 3);
  BOOST_CHECK_EQUAL(top.size(), 3);
  BOOST_CHECK_EQUAL(bottom[0], 0);
  BOOST_CHECK_EQUAL(bottom[1], 1);
  BOOST_CHECK_EQUAL(bottom[2], 2);
  BOOST_CHECK_EQUAL(top[0], 4);
  BOOST_CHECK_EQUAL(top[1], 5);
  BOOST_CHECK_EQUAL(top[2], 6);
  BOOST_CHECK(b3.get_bounding_corners().first == bottom);
  BOOST_CHECK(b3.get_bounding_corners().second == top);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(box_other, T, list_of_tested_variants)
{
  using P = typename Box<T>::Point_t;
  P nan = {
      std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};
  P inf = {P::T_inf, P::T_inf, P::T_inf};
  P minus_inf = {-P::T_inf, -P::T_inf, -P::T_inf};

  Box<T> b({0, 1, 2}, {4, 5, 6});

  BOOST_CHECK(b.contains({1, 2, 3}));
  BOOST_CHECK(!b.contains({1, 8, 3}));
  if (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(!b.contains(nan));
  BOOST_CHECK(!b.contains(inf));
  BOOST_CHECK(!b.contains(minus_inf));

  BOOST_CHECK_EQUAL(b.dimension(), 3);

  b.inflate(2);
  auto& bottom = b.get_lower_corner();
  auto& top = b.get_upper_corner();
  BOOST_CHECK_EQUAL(bottom[0], -2);
  BOOST_CHECK_EQUAL(bottom[1], -1);
  BOOST_CHECK_EQUAL(bottom[2], 0);
  BOOST_CHECK_EQUAL(top[0], 6);
  BOOST_CHECK_EQUAL(top[1], 7);
  BOOST_CHECK_EQUAL(top[2], 8);

  top = inf;

  BOOST_CHECK(b.contains({1, 2, 3}));
  BOOST_CHECK(b.contains({1, 8, 3}));
  if (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(!b.contains(nan));
  BOOST_CHECK(b.contains(inf));
  BOOST_CHECK(!b.contains(minus_inf));

  BOOST_CHECK_EQUAL(b.dimension(), 3);

  b.inflate(2);
  BOOST_CHECK_EQUAL(bottom[0], -4);
  BOOST_CHECK_EQUAL(bottom[1], -3);
  BOOST_CHECK_EQUAL(bottom[2], -2);
  BOOST_CHECK_EQUAL(top, inf);

  bottom = minus_inf;
  top = {4, 5, 6};

  BOOST_CHECK(b.contains({1, 2, 3}));
  BOOST_CHECK(!b.contains({1, 8, 3}));
  if (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(!b.contains(nan));
  BOOST_CHECK(!b.contains(inf));
  BOOST_CHECK(b.contains(minus_inf));

  BOOST_CHECK_EQUAL(b.dimension(), 3);

  b.inflate(2);
  BOOST_CHECK_EQUAL(bottom, minus_inf);
  BOOST_CHECK_EQUAL(top[0], 6);
  BOOST_CHECK_EQUAL(top[1], 7);
  BOOST_CHECK_EQUAL(top[2], 8);

  top = inf;

  BOOST_CHECK(b.contains({1, 2, 3}));
  BOOST_CHECK(b.contains({1, 8, 3}));
  if (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(!b.contains(nan));
  BOOST_CHECK(b.contains(inf));
  BOOST_CHECK(b.contains(minus_inf));

  BOOST_CHECK_EQUAL(b.dimension(), 3);

  b.inflate(2);
  BOOST_CHECK_EQUAL(bottom, minus_inf);
  BOOST_CHECK_EQUAL(top, inf);
}
