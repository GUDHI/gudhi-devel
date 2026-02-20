/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
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

using list_of_tested_variants = boost::mpl::list<double, float, int>;

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
  P minus_inf = {P::T_m_inf, P::T_m_inf, P::T_m_inf};

  Box<T> b({0, 1, 2}, {4, 5, 6});

  BOOST_CHECK(b.contains({1, 2, 3}));
  BOOST_CHECK(!b.contains({1, 8, 3}));
  if (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(!b.contains(nan));
  BOOST_CHECK(!b.contains(inf));
  BOOST_CHECK(!b.contains(minus_inf));

  BOOST_CHECK_EQUAL(b.get_dimension(), 3);

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

  BOOST_CHECK_EQUAL(b.get_dimension(), 3);

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

  BOOST_CHECK_EQUAL(b.get_dimension(), 3);

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

  BOOST_CHECK_EQUAL(b.get_dimension(), 3);

  b.inflate(2);
  BOOST_CHECK_EQUAL(bottom, minus_inf);
  BOOST_CHECK_EQUAL(top, inf);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(box_friends, T, list_of_tested_variants)
{
  using P = typename Box<T>::Point_t;
  P inf = {P::T_inf, P::T_inf, P::T_inf};
  P minus_inf = {P::T_m_inf, P::T_m_inf, P::T_m_inf};

  Box<T> a({0, 1, 2}, {4, 5, 6});
  Box<T> b({2, 1, 0}, {6, 5, 4});

  auto res = get_smallest_enclosing_box(a, b);
  auto bottom = res.get_lower_corner();
  auto top = res.get_upper_corner();
  BOOST_CHECK_EQUAL(bottom.size(), 3);
  BOOST_CHECK_EQUAL(top.size(), 3);
  BOOST_CHECK_EQUAL(bottom[0], 0);
  BOOST_CHECK_EQUAL(bottom[1], 1);
  BOOST_CHECK_EQUAL(bottom[2], 0);
  BOOST_CHECK_EQUAL(top[0], 6);
  BOOST_CHECK_EQUAL(top[1], 5);
  BOOST_CHECK_EQUAL(top[2], 6);

  res = get_smallest_enclosing_box(Box(inf, inf), Box(minus_inf, minus_inf));
  bottom = res.get_lower_corner();
  top = res.get_upper_corner();
  BOOST_CHECK_EQUAL(bottom.size(), 0);
  BOOST_CHECK_EQUAL(top.size(), 0);

  res = get_smallest_enclosing_box(a, Box(inf, inf));
  bottom = res.get_lower_corner();
  top = res.get_upper_corner();
  BOOST_CHECK_EQUAL(bottom.size(), 3);
  BOOST_CHECK_EQUAL(top.size(), 3);
  BOOST_CHECK_EQUAL(bottom[0], 0);
  BOOST_CHECK_EQUAL(bottom[1], 1);
  BOOST_CHECK_EQUAL(bottom[2], 2);
  BOOST_CHECK_EQUAL(top[0], 4);
  BOOST_CHECK_EQUAL(top[1], 5);
  BOOST_CHECK_EQUAL(top[2], 6);

  res = get_smallest_enclosing_box(a, Box({4, 5, 6}, inf));
  bottom = res.get_lower_corner();
  top = res.get_upper_corner();
  BOOST_CHECK_EQUAL(bottom.size(), 3);
  BOOST_CHECK_EQUAL(top.size(), 3);
  BOOST_CHECK_EQUAL(bottom[0], 0);
  BOOST_CHECK_EQUAL(bottom[1], 1);
  BOOST_CHECK_EQUAL(bottom[2], 2);
  BOOST_CHECK_EQUAL(top[0], P::T_inf);
  BOOST_CHECK_EQUAL(top[1], P::T_inf);
  BOOST_CHECK_EQUAL(top[2], P::T_inf);

  res = get_smallest_enclosing_box(a, Box(minus_inf, {0, 1, 2}));
  bottom = res.get_lower_corner();
  top = res.get_upper_corner();
  BOOST_CHECK_EQUAL(bottom.size(), 3);
  BOOST_CHECK_EQUAL(top.size(), 3);
  BOOST_CHECK_EQUAL(bottom[0], P::T_m_inf);
  BOOST_CHECK_EQUAL(bottom[1], P::T_m_inf);
  BOOST_CHECK_EQUAL(bottom[2], P::T_m_inf);
  BOOST_CHECK_EQUAL(top[0], 4);
  BOOST_CHECK_EQUAL(top[1], 5);
  BOOST_CHECK_EQUAL(top[2], 6);
}
