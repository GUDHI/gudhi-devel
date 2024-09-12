/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/Line.h>

using Gudhi::multi_persistence::Line;

typedef boost::mpl::list<double, float, int> list_of_tested_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(line_constructors, T, list_of_tested_variants)
{
  using P = typename Line<T>::Point;

  Line<T> l;
  BOOST_CHECK(l.base_point().empty());
  BOOST_CHECK(l.direction().empty());

  P x({1,2,3});

  Line<T> l2(x);
  auto bp = l2.base_point();
  BOOST_CHECK(bp.size() == 3);
  BOOST_CHECK_EQUAL(bp[0], 1);
  BOOST_CHECK_EQUAL(bp[1], 2);
  BOOST_CHECK_EQUAL(bp[2], 3);
  BOOST_CHECK(l2.direction().empty());

  Line<T> l3(std::move(x));
  bp = l3.base_point();
  BOOST_CHECK(bp.size() == 3);
  BOOST_CHECK_EQUAL(bp[0], 1);
  BOOST_CHECK_EQUAL(bp[1], 2);
  BOOST_CHECK_EQUAL(bp[2], 3);
  BOOST_CHECK(l3.direction().empty());

  Line<T> l4({1,2,3}, {4,5,6});
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

  BOOST_CHECK_THROW(Line<T> l5({1,2,3}, {4,-5,6}), std::invalid_argument);
  BOOST_CHECK_THROW(Line<T> l6({1,2,3}, {0,0,0}), std::invalid_argument);
  BOOST_CHECK_THROW(Line<T> l7({1,2,3}, {1,2,3,4}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(line_intersections, T, list_of_tested_variants)
{
  using P = typename Line<T>::Point;
  using KP = typename Line<T>::K_critical_point;

  const double inf = std::numeric_limits<double>::infinity();

  Line<T> l({1,2,3}, {4,0,6});

  double t = l.template compute_forward_intersection<double>(P{2,1,3});
  BOOST_CHECK_EQUAL(t, 0.25);
  t = l.template compute_forward_intersection<double>(P{2,3,3});
  BOOST_CHECK_EQUAL(t, inf);

  t = l.template compute_forward_intersection<double>(KP({P{2,1,3}, P{2,3,3}}));
  BOOST_CHECK_EQUAL(t, 0.25);

  t = l.template compute_backward_intersection<double>(P{2,3,3});
  BOOST_CHECK_EQUAL(t, 0);
  t = l.template compute_backward_intersection<double>(P{2,1,3});
  BOOST_CHECK_EQUAL(t, -inf);

  t = l.template compute_backward_intersection<double>(KP({P{2,1,3}, P{2,3,3}}));
  BOOST_CHECK_EQUAL(t, 0);

  std::pair<P, P> bounds = l.get_bounds({{-10, 0, 10}, {10, 4, 10}});
  auto& bottom = bounds.first;
  auto& top = bounds.second;
  BOOST_CHECK_EQUAL(bottom.size(), 3);
  BOOST_CHECK_EQUAL(bottom[0], T(5. + 2. / 3.));
  BOOST_CHECK_EQUAL(bottom[1], 2);
  BOOST_CHECK_EQUAL(bottom[2], T(3. + T(7. / 6.) * 6.));
  BOOST_CHECK_EQUAL(top.size(), 3);
  BOOST_CHECK_EQUAL(top[0], T(5. + 2. / 3.));
  BOOST_CHECK_EQUAL(top[1], 2);
  BOOST_CHECK_EQUAL(top[2], T(3. + T(7. / 6.) * 6.));

  bounds = l.get_bounds({{-10, 0, 10}, {10, 1, 10}});
  BOOST_CHECK(bounds.first.is_plus_inf());
  BOOST_CHECK(bounds.second.is_plus_inf());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(line_other, T, list_of_tested_variants)
{
  using P = typename Line<T>::Point;

  Line<T> l({1,2,3}, {4,0,6});

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



