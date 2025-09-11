/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024-25 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/Point.h>
#include <gudhi/Multi_persistence/Line.h>
#include <gudhi/Multi_parameter_filtration.h>

using Gudhi::multi_filtration::Multi_parameter_filtration;
using Gudhi::multi_persistence::Line;
using Gudhi::multi_persistence::Point;

typedef boost::mpl::list<double, float, int> list_of_tested_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(line_constructors, T, list_of_tested_variants)
{
  using P = typename Line<T>::Point_t;

  Line<T> l;
  BOOST_CHECK(l.base_point().size() == 0);
  BOOST_CHECK(l.direction().size() == 0);

  P x({1, 2, 3});

  Line<T> l2(x);
  auto bp = l2.base_point();
  BOOST_CHECK(bp.size() == 3);
  BOOST_CHECK_EQUAL(bp[0], 1);
  BOOST_CHECK_EQUAL(bp[1], 2);
  BOOST_CHECK_EQUAL(bp[2], 3);
  BOOST_CHECK(l2.direction().size() == 0);

  Line<T> l3(std::move(x));
  bp = l3.base_point();
  BOOST_CHECK(bp.size() == 3);
  BOOST_CHECK_EQUAL(bp[0], 1);
  BOOST_CHECK_EQUAL(bp[1], 2);
  BOOST_CHECK_EQUAL(bp[2], 3);
  BOOST_CHECK(l3.direction().size() == 0);

  Line<T> l4({1, 2, 3}, {4, 5, 6});
  bp = l4.base_point();
  auto dir = l4.direction();
  BOOST_CHECK(bp.size() == 3);
  BOOST_CHECK_EQUAL(bp[0], 1);
  BOOST_CHECK_EQUAL(bp[1], 2);
  BOOST_CHECK_EQUAL(bp[2], 3);
  BOOST_CHECK(dir.size() == 3);
  BOOST_CHECK_EQUAL(dir[0], 4);
  BOOST_CHECK_EQUAL(dir[1], 5);
  BOOST_CHECK_EQUAL(dir[2], 6);

  BOOST_CHECK_THROW(Line<T> l5({1, 2, 3}, {4, -5, 6}), std::invalid_argument);
  BOOST_CHECK_THROW(Line<T> l6({1, 2, 3}, {0, 0, 0}), std::invalid_argument);
  BOOST_CHECK_THROW(Line<T> l7({1, 2, 3}, {1, 2, 3, 4}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(line_intersections, T, list_of_tested_variants)
{
  using F = Multi_parameter_filtration<T>;
  Point<T> p_inf = {Point<T>::T_inf, Point<T>::T_inf, Point<T>::T_inf};
  Point<T> p_minus_inf = {Point<T>::T_m_inf, Point<T>::T_m_inf, Point<T>::T_m_inf};

  const double inf = std::numeric_limits<double>::infinity();

  Line<T> l({1, 2, 3}, {4, 0, 6});

  double t = l.template compute_forward_intersection<double>(F({2, 1, 3}));
  BOOST_CHECK_EQUAL(t, 0.25);
  t = l.template compute_forward_intersection<double>(F({2, 3, 3}));
  BOOST_CHECK_EQUAL(t, inf);

  t = l.template compute_forward_intersection<double>(F({2, 1, 3, 2, 3, 3}, 3));
  BOOST_CHECK_EQUAL(t, 0.25);

  t = l.template compute_backward_intersection<double>(F({2, 3, 3}));
  BOOST_CHECK_EQUAL(t, 0);
  t = l.template compute_backward_intersection<double>(F({2, 1, 3}));
  BOOST_CHECK_EQUAL(t, -inf);

  t = l.template compute_backward_intersection<double>(F({2, 1, 3, 2, 3, 3}, 3));
  BOOST_CHECK_EQUAL(t, 0);

  std::pair<T, T> bounds = l.get_bounds({{-10, 0, 10}, {10, 4, 10}});
  auto bottom = l[bounds.first];
  auto top = l[bounds.second];
  BOOST_CHECK_EQUAL(bottom.size(), 3);
  BOOST_CHECK_EQUAL(bottom[0], T(5. + 2. / 3.));
  BOOST_CHECK_EQUAL(bottom[1], 2);
  BOOST_CHECK_EQUAL(bottom[2], T(3. + T(7. / 6.) * 6.));
  BOOST_CHECK_EQUAL(top.size(), 3);
  BOOST_CHECK_EQUAL(top[0], T(5. + 2. / 3.));
  BOOST_CHECK_EQUAL(top[1], 2);
  BOOST_CHECK_EQUAL(top[2], T(3. + T(7. / 6.) * 6.));

  bounds = l.get_bounds({{-10, 0, 10}, {10, 1, 10}});
  BOOST_CHECK_EQUAL(bounds.first, F::T_inf);
  BOOST_CHECK_EQUAL(bounds.second, F::T_m_inf);
  BOOST_CHECK(l[bounds.first] == p_inf);
  BOOST_CHECK(l[bounds.second] == p_minus_inf);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(line_other, T, list_of_tested_variants)
{
  using P = typename Line<T>::Point_t;

  Line<T> l({1, 2, 3}, {4, 0, 6});

  P x = l[1.25];
  BOOST_CHECK_EQUAL(x.size(), 3);
  BOOST_CHECK_EQUAL(x[0], 1. + T(1.25) * 4.);
  BOOST_CHECK_EQUAL(x[1], 2);
  BOOST_CHECK_EQUAL(x[2], 3. + T(1.25) * 6.);

  x = l[0];
  BOOST_CHECK_EQUAL(x.size(), 3);
  BOOST_CHECK_EQUAL(x[0], 1);
  BOOST_CHECK_EQUAL(x[1], 2);
  BOOST_CHECK_EQUAL(x[2], 3);

  x = l[-3.25];
  BOOST_CHECK_EQUAL(x.size(), 3);
  BOOST_CHECK_EQUAL(x[0], 1. + T(-3.25) * 4.);
  BOOST_CHECK_EQUAL(x[1], 2);
  BOOST_CHECK_EQUAL(x[2], 3. + T(-3.25) * 6.);

  l += {2, 5, 6};
  x = l[1.25];
  BOOST_CHECK_EQUAL(x.size(), 3);
  BOOST_CHECK_EQUAL(x[0], 3. + T(1.25) * 4.);
  BOOST_CHECK_EQUAL(x[1], 7);
  BOOST_CHECK_EQUAL(x[2], 9. + T(1.25) * 6.);

  x = l[0];
  BOOST_CHECK_EQUAL(x.size(), 3);
  BOOST_CHECK_EQUAL(x[0], 3);
  BOOST_CHECK_EQUAL(x[1], 7);
  BOOST_CHECK_EQUAL(x[2], 9);

  x = l[-3.25];
  BOOST_CHECK_EQUAL(x.size(), 3);
  BOOST_CHECK_EQUAL(x[0], 3. + T(-3.25) * 4.);
  BOOST_CHECK_EQUAL(x[1], 7);
  BOOST_CHECK_EQUAL(x[2], 9. + T(-3.25) * 6.);
}
