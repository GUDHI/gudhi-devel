/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <boost/test/tools/old/interface.hpp>
#include <cmath>
#include <limits>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_filtration"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/One_critical_filtration.h>

using Gudhi::multi_filtration::One_critical_filtration;

typedef boost::mpl::list<double, float, int> list_of_tested_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(one_critical_filtration_constructors, T, list_of_tested_variants)
{
  One_critical_filtration<T> f;
  BOOST_CHECK(!f.empty());
  BOOST_CHECK(f.num_parameters() == 1);

  One_critical_filtration<T> f1(3);
  BOOST_CHECK(f1.size() == 3);
  BOOST_CHECK(f1.num_parameters() == 3);
  BOOST_CHECK(f1[0] == -One_critical_filtration<T>::T_inf);

  One_critical_filtration<T> f2(3, 0);
  BOOST_CHECK(f2.size() == 3);
  BOOST_CHECK(f2.num_parameters() == 3);
  BOOST_CHECK(f2[0] == 0);

  One_critical_filtration<T> f3({0, 1, 2});
  BOOST_CHECK(f3.size() == 3);
  BOOST_CHECK(f3.num_parameters() == 3);
  BOOST_CHECK(f3[0] == 0);
  BOOST_CHECK(f3[1] == 1);
  BOOST_CHECK(f3[2] == 2);

  std::vector<T> v{0, 1, 2};
  One_critical_filtration<T> f4(v);
  BOOST_CHECK(f4.size() == 3);
  BOOST_CHECK(f4.num_parameters() == 3);
  BOOST_CHECK(f4[0] == 0);
  BOOST_CHECK(f4[1] == 1);
  BOOST_CHECK(f4[2] == 2);

  One_critical_filtration<T> f5(v.begin(), v.end());
  BOOST_CHECK(f5.size() == 3);
  BOOST_CHECK(f5.num_parameters() == 3);
  BOOST_CHECK(f5[0] == 0);
  BOOST_CHECK(f5[1] == 1);
  BOOST_CHECK(f5[2] == 2);

  f = f5;
  BOOST_CHECK(f.size() == 3);
  BOOST_CHECK(f.num_parameters() == 3);
  BOOST_CHECK(f[0] == 0);
  BOOST_CHECK(f[1] == 1);
  BOOST_CHECK(f[2] == 2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_critical_filtration_utilities, T, list_of_tested_variants)
{
  One_critical_filtration<T> f({0, 1, 2});
  bool test = std::is_same_v<decltype(f[0]), T&>;
  BOOST_CHECK(test);

  One_critical_filtration<float> f2 = f.template as_type<float>();
  test = std::is_same_v<decltype(f2[0]), float&>;
  BOOST_CHECK(test);
  BOOST_CHECK(f2.size() == 3);
  BOOST_CHECK(f2.num_parameters() == 3);
  BOOST_CHECK(f2[0] == 0.);
  BOOST_CHECK(f2[1] == 1.);
  BOOST_CHECK(f2[2] == 2.);

  BOOST_CHECK(!f.is_inf());
  BOOST_CHECK(!f.is_minus_inf());
  BOOST_CHECK(!f.is_nan());
  BOOST_CHECK(f.is_finite());

  One_critical_filtration<T> f3;
  BOOST_CHECK(!f3.is_inf());
  BOOST_CHECK(f3.is_minus_inf());
  BOOST_CHECK(!f3.is_nan());
  BOOST_CHECK(!f3.is_finite());

  //{-inf, -inf, -inf} is considered finite as the user is supposed to updates the values to something else
  //the idea is just to reserve space and to give the possibility to use `f4[i] =`
  //if the value should really be -inf, use `f4(1)` or `f4 = minus_inf()` instead.
  One_critical_filtration<T> f4(3);
  BOOST_CHECK(!f4.is_inf());
  BOOST_CHECK(!f4.is_minus_inf());
  BOOST_CHECK(!f4.is_nan());
  BOOST_CHECK(f4.is_finite());

  One_critical_filtration<T> f5(1);
  BOOST_CHECK(!f5.is_inf());
  BOOST_CHECK(f5.is_minus_inf());
  BOOST_CHECK(!f5.is_nan());
  BOOST_CHECK(!f5.is_finite());

  auto f6 = One_critical_filtration<T>::nan();
  BOOST_CHECK(!f6.is_inf());
  BOOST_CHECK(!f6.is_minus_inf());
  BOOST_CHECK(f6.is_nan());
  BOOST_CHECK(!f6.is_finite());

  auto f7 = One_critical_filtration<T>::inf();
  BOOST_CHECK(f7.is_inf());
  BOOST_CHECK(!f7.is_minus_inf());
  BOOST_CHECK(!f7.is_nan());
  BOOST_CHECK(!f7.is_finite());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_critical_filtration_comparators, T, list_of_tested_variants)
{
  One_critical_filtration<T> f1({0, 1, 2});
  One_critical_filtration<T> f2({-1, 0, 1});
  One_critical_filtration<T> f3({0, 2, 3});
  One_critical_filtration<T> f4({5, -1, 2});

  BOOST_CHECK(!(f1 < f1));
  BOOST_CHECK(!(f1 < f2));
  BOOST_CHECK(f1 < f3);
  BOOST_CHECK(!(f1 < f4));
  BOOST_CHECK(!(f1 < One_critical_filtration<T>::nan()));
  BOOST_CHECK(f1 < One_critical_filtration<T>::inf());
  BOOST_CHECK(!(f1 < One_critical_filtration<T>::minus_inf()));

  BOOST_CHECK(f1 <= f1);
  BOOST_CHECK(!(f1 <= f2));
  BOOST_CHECK(f1 <= f3);
  BOOST_CHECK(!(f1 <= f4));
  BOOST_CHECK(!(f1 <= One_critical_filtration<T>::nan()));
  BOOST_CHECK(f1 <= One_critical_filtration<T>::inf());
  BOOST_CHECK(!(f1 <= One_critical_filtration<T>::minus_inf()));

  BOOST_CHECK(!(f1 > f1));
  BOOST_CHECK(f1 > f2);
  BOOST_CHECK(!(f1 > f3));
  BOOST_CHECK(!(f1 > f4));
  BOOST_CHECK(!(f1 > One_critical_filtration<T>::nan()));
  BOOST_CHECK(!(f1 > One_critical_filtration<T>::inf()));
  BOOST_CHECK(f1 > One_critical_filtration<T>::minus_inf());

  BOOST_CHECK(f1 >= f1);
  BOOST_CHECK(f1 >= f2);
  BOOST_CHECK(!(f1 >= f3));
  BOOST_CHECK(!(f1 >= f4));
  BOOST_CHECK(!(f1 >= One_critical_filtration<T>::nan()));
  BOOST_CHECK(!(f1 >= One_critical_filtration<T>::inf()));
  BOOST_CHECK(f1 >= One_critical_filtration<T>::minus_inf());

  BOOST_CHECK(f1 == f1);
  BOOST_CHECK(!(f1 == f2));
  BOOST_CHECK(!(f1 == f3));
  BOOST_CHECK(!(f1 == f4));
  BOOST_CHECK(!(f1 == One_critical_filtration<T>::nan()));
  BOOST_CHECK(!(f1 == One_critical_filtration<T>::inf()));
  BOOST_CHECK(!(f1 == One_critical_filtration<T>::minus_inf()));

  BOOST_CHECK(!(f1 != f1));
  BOOST_CHECK(f1 != f2);
  BOOST_CHECK(f1 != f3);
  BOOST_CHECK(f1 != f4);
  BOOST_CHECK(f1 != One_critical_filtration<T>::nan());
  BOOST_CHECK(f1 != One_critical_filtration<T>::inf());
  BOOST_CHECK(f1 != One_critical_filtration<T>::minus_inf());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_critical_filtration_operators, T, list_of_tested_variants)
{;
  One_critical_filtration<T> f({-10, 0, 1});
  One_critical_filtration<T> f2({5, 2, -1});
  One_critical_filtration<T> f3(
      {-One_critical_filtration<T>::T_inf, One_critical_filtration<T>::T_inf, -One_critical_filtration<T>::T_inf});
  One_critical_filtration<T> f4(
      {One_critical_filtration<T>::T_inf, -One_critical_filtration<T>::T_inf, One_critical_filtration<T>::T_inf});

  auto res = -f;
  BOOST_CHECK_EQUAL(res[0], 10);
  BOOST_CHECK_EQUAL(res[1], 0);
  BOOST_CHECK_EQUAL(res[2], -1);
  BOOST_CHECK((-One_critical_filtration<T>::inf()).is_minus_inf());
  BOOST_CHECK((-One_critical_filtration<T>::minus_inf()).is_inf());
  BOOST_CHECK((-One_critical_filtration<T>::nan()).is_nan());

  res = f - f2;
  BOOST_CHECK_EQUAL(res[0], -15);
  BOOST_CHECK_EQUAL(res[1], -2);
  BOOST_CHECK_EQUAL(res[2], 2);

  res = f - f3;
  BOOST_CHECK_EQUAL(res[0], f4[0]);
  BOOST_CHECK_EQUAL(res[1], f4[1]);
  BOOST_CHECK_EQUAL(res[2], f4[2]);

  res = f3 - f;
  BOOST_CHECK_EQUAL(res[0], f3[0]);
  BOOST_CHECK_EQUAL(res[1], f3[1]);
  BOOST_CHECK_EQUAL(res[2], f3[2]);

  res = 5 - f;
  BOOST_CHECK_EQUAL(res[0], 15);
  BOOST_CHECK_EQUAL(res[1], 5);
  BOOST_CHECK_EQUAL(res[2], 4);

  res = f - 5;
  BOOST_CHECK_EQUAL(res[0], -15);
  BOOST_CHECK_EQUAL(res[1], -5);
  BOOST_CHECK_EQUAL(res[2], -4);

  BOOST_CHECK((f - One_critical_filtration<T>::inf()).is_minus_inf());
  BOOST_CHECK((One_critical_filtration<T>::inf() - f).is_inf());
  BOOST_CHECK((f - One_critical_filtration<T>::minus_inf()).is_inf());
  BOOST_CHECK((One_critical_filtration<T>::minus_inf() - f).is_minus_inf());
  BOOST_CHECK((f - One_critical_filtration<T>::nan()).is_nan());
  BOOST_CHECK((One_critical_filtration<T>::nan() - f).is_nan());

  res = f3 - f3;
  BOOST_CHECK(res.is_nan());
  res = f3 - f4;
  BOOST_CHECK_EQUAL(res[0], f3[0]);
  BOOST_CHECK_EQUAL(res[1], f3[1]);
  BOOST_CHECK_EQUAL(res[2], f3[2]);

  res = f + f2;
  BOOST_CHECK_EQUAL(res[0], -5);
  BOOST_CHECK_EQUAL(res[1], 2);
  BOOST_CHECK_EQUAL(res[2], 0);

  res = f + f3;
  BOOST_CHECK_EQUAL(res[0], f3[0]);
  BOOST_CHECK_EQUAL(res[1], f3[1]);
  BOOST_CHECK_EQUAL(res[2], f3[2]);

  res = 5 + f;
  BOOST_CHECK_EQUAL(res[0], -5);
  BOOST_CHECK_EQUAL(res[1], 5);
  BOOST_CHECK_EQUAL(res[2], 6);

  res = f + 5;
  BOOST_CHECK_EQUAL(res[0], -5);
  BOOST_CHECK_EQUAL(res[1], 5);
  BOOST_CHECK_EQUAL(res[2], 6);

  BOOST_CHECK((f + One_critical_filtration<T>::inf()).is_inf());
  BOOST_CHECK((One_critical_filtration<T>::inf() + f).is_inf());
  BOOST_CHECK((f + One_critical_filtration<T>::minus_inf()).is_minus_inf());
  BOOST_CHECK((One_critical_filtration<T>::minus_inf() + f).is_minus_inf());
  BOOST_CHECK((f + One_critical_filtration<T>::nan()).is_nan());
  BOOST_CHECK((One_critical_filtration<T>::nan() + f).is_nan());

  res = f3 + f4;
  BOOST_CHECK(res.is_nan());
  res = f3 + f3;
  BOOST_CHECK_EQUAL(res[0], f3[0]);
  BOOST_CHECK_EQUAL(res[1], f3[1]);
  BOOST_CHECK_EQUAL(res[2], f3[2]);

  res = f * f2;
  BOOST_CHECK_EQUAL(res[0], -50);
  BOOST_CHECK_EQUAL(res[1], 0);
  BOOST_CHECK_EQUAL(res[2], -1);

  res = f * f3;
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(res[0], f4[0]);
    BOOST_CHECK(std::isnan(res[1]));
    BOOST_CHECK_EQUAL(res[2], f3[2]);
  } else {
    BOOST_CHECK(res.is_nan());
  }

  res = 5 * f;
  BOOST_CHECK_EQUAL(res[0], -50);
  BOOST_CHECK_EQUAL(res[1], 0);
  BOOST_CHECK_EQUAL(res[2], 5);

  res = f * 5;
  BOOST_CHECK_EQUAL(res[0], -50);
  BOOST_CHECK_EQUAL(res[1],  0);
  BOOST_CHECK_EQUAL(res[2], 5);

  res = f * One_critical_filtration<T>::inf();
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(res[0], -One_critical_filtration<T>::T_inf);
    BOOST_CHECK(std::isnan(res[1]));
    BOOST_CHECK_EQUAL(res[2], One_critical_filtration<T>::T_inf);
  } else {
    BOOST_CHECK(res.is_nan());
  }
  res = One_critical_filtration<T>::inf() * f;
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(res[0], -One_critical_filtration<T>::T_inf);
    BOOST_CHECK(std::isnan(res[1]));
    BOOST_CHECK_EQUAL(res[2], One_critical_filtration<T>::T_inf);
  } else {
    BOOST_CHECK(res.is_nan());
  }
  res = f * One_critical_filtration<T>::minus_inf();
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(res[0], One_critical_filtration<T>::T_inf);
    BOOST_CHECK(std::isnan(res[1]));
    BOOST_CHECK_EQUAL(res[2], -One_critical_filtration<T>::T_inf);
  } else {
    BOOST_CHECK(res.is_nan());
  }
  res = One_critical_filtration<T>::minus_inf() * f;
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(res[0], One_critical_filtration<T>::T_inf);
    BOOST_CHECK(std::isnan(res[1]));
    BOOST_CHECK_EQUAL(res[2], -One_critical_filtration<T>::T_inf);
  } else {
    BOOST_CHECK(res.is_nan());
  }
  res = f * One_critical_filtration<T>::nan();
  BOOST_CHECK(res.is_nan());
  res = One_critical_filtration<T>::nan() * f;
  BOOST_CHECK(res.is_nan());

  res = f3 * f3;
  BOOST_CHECK(res.is_inf());
  res = f3 * f4;
  BOOST_CHECK(res.is_minus_inf());

  res = f / f2;
  BOOST_CHECK_EQUAL(res[0], -2);
  BOOST_CHECK_EQUAL(res[1], 0);
  BOOST_CHECK_EQUAL(res[2], -1);

  res = f / f3;
  BOOST_CHECK_EQUAL(res[0], 0);
  BOOST_CHECK_EQUAL(res[1], 0);
  BOOST_CHECK_EQUAL(res[2], 0);

  res = f3 / f;
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(res[0], f4[0]);
    BOOST_CHECK(std::isnan(res[1]));
    BOOST_CHECK_EQUAL(res[2], f3[2]);
  } else {
    BOOST_CHECK(res.is_nan());
  }

  res = 5 / f;
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(res[0], -0.5);
    BOOST_CHECK(std::isnan(res[1]));
    BOOST_CHECK_EQUAL(res[2], 5);
  } else {
    BOOST_CHECK(res.is_nan());
  }

  res = f / 5;
  BOOST_CHECK_EQUAL(res[0], -2);
  BOOST_CHECK_EQUAL(res[1],  0);
  BOOST_CHECK_EQUAL(res[2], static_cast<T>(1) / static_cast<T>(5)); //to avoid precision error

  res = f / One_critical_filtration<T>::inf();
  BOOST_CHECK_EQUAL(res[0], 0);
  BOOST_CHECK_EQUAL(res[1], 0);
  BOOST_CHECK_EQUAL(res[2], 0);
  res = One_critical_filtration<T>::inf() / f;
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(res[0], -One_critical_filtration<T>::T_inf);
    BOOST_CHECK(std::isnan(res[1]));
    BOOST_CHECK_EQUAL(res[2], One_critical_filtration<T>::T_inf);
  } else {
    BOOST_CHECK(res.is_nan());
  }
  res = f / One_critical_filtration<T>::minus_inf();
  BOOST_CHECK_EQUAL(res[0], 0);
  BOOST_CHECK_EQUAL(res[1], 0);
  BOOST_CHECK_EQUAL(res[2], 0);
  res = One_critical_filtration<T>::minus_inf() / f;
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(res[0], One_critical_filtration<T>::T_inf);
    BOOST_CHECK(std::isnan(res[1]));
    BOOST_CHECK_EQUAL(res[2], -One_critical_filtration<T>::T_inf);
  } else {
    BOOST_CHECK(res.is_nan());
  }
  res = f / One_critical_filtration<T>::nan();
  BOOST_CHECK(res.is_nan());
  res = One_critical_filtration<T>::nan() / f;
  BOOST_CHECK(res.is_nan());

  res = f3 / f3;
  BOOST_CHECK(res.is_nan());
  res = f3 / f4;
  BOOST_CHECK(res.is_nan());
  res = f / One_critical_filtration<T>({0, 0, 0});
  BOOST_CHECK(res.is_nan());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_critical_filtration_modifiers, T, list_of_tested_variants)
{
  One_critical_filtration<T> f({0, 1, 2});
  BOOST_CHECK_EQUAL(f[0], 0);
  BOOST_CHECK_EQUAL(f[1], 1);
  BOOST_CHECK_EQUAL(f[2], 2);

  f.push_to_least_common_upper_bound({-1, 5, 6});
  BOOST_CHECK_EQUAL(f[0], 0);
  BOOST_CHECK_EQUAL(f[1], 5);
  BOOST_CHECK_EQUAL(f[2], 6);

  f.push_to_least_common_upper_bound({-1, -5, -6});
  BOOST_CHECK_EQUAL(f[0], 0);
  BOOST_CHECK_EQUAL(f[1], 5);
  BOOST_CHECK_EQUAL(f[2], 6);

  f.push_to_least_common_upper_bound(One_critical_filtration<T>::minus_inf());
  BOOST_CHECK_EQUAL(f[0], 0);
  BOOST_CHECK_EQUAL(f[1], 5);
  BOOST_CHECK_EQUAL(f[2], 6);

  f.push_to_least_common_upper_bound(One_critical_filtration<T>::inf());
  BOOST_CHECK(f.is_inf());

  f.push_to_least_common_upper_bound(One_critical_filtration<T>::nan());
  BOOST_CHECK(f.is_inf());

  f.pull_to_greatest_common_lower_bound({-1, 5, 6});
  BOOST_CHECK_EQUAL(f[0], -1);
  BOOST_CHECK_EQUAL(f[1], 5);
  BOOST_CHECK_EQUAL(f[2], 6);

  f.pull_to_greatest_common_lower_bound({1, 8, 9});
  BOOST_CHECK_EQUAL(f[0], -1);
  BOOST_CHECK_EQUAL(f[1], 5);
  BOOST_CHECK_EQUAL(f[2], 6);

  f.pull_to_greatest_common_lower_bound(One_critical_filtration<T>::inf());
  BOOST_CHECK_EQUAL(f[0], -1);
  BOOST_CHECK_EQUAL(f[1], 5);
  BOOST_CHECK_EQUAL(f[2], 6);

  f.pull_to_greatest_common_lower_bound(One_critical_filtration<T>::minus_inf());
  BOOST_CHECK(f.is_minus_inf());

  f.pull_to_greatest_common_lower_bound(One_critical_filtration<T>::nan());
  BOOST_CHECK(f.is_minus_inf());

  std::vector<std::vector<int> > grid = {{0, 2, 4, 8}, {0, 3, 6, 9}, {0, 4, 8, 16}};

  f.push_to_least_common_upper_bound({1, 7, 5});
  f.project_onto_grid(grid, true);
  BOOST_CHECK_EQUAL(f[0], 1);
  BOOST_CHECK_EQUAL(f[1], 3);
  BOOST_CHECK_EQUAL(f[2], 2);

  f.push_to_least_common_upper_bound({1, 7, 5});
  f.project_onto_grid(grid, false);
  BOOST_CHECK_EQUAL(f[0], 2);
  BOOST_CHECK_EQUAL(f[1], 9);
  BOOST_CHECK_EQUAL(f[2], 8);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_critical_filtration_friends, T, list_of_tested_variants)
{
  One_critical_filtration<T> f({0, 1, 2});
  
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {2,3,5,9}), 13);
  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(5))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, {2,3,5}), static_cast<T>(std::sqrt(T(17))));

  f = {1, 7, 5};

  std::vector<std::vector<int> > grid = {{0, 2, 4, 8}, {0, 3, 6, 9}, {0, 4, 8, 16}};
  auto res = compute_coordinates_in_grid(f, grid);
  BOOST_CHECK_EQUAL(res[0], 1);
  BOOST_CHECK_EQUAL(res[1], 3);
  BOOST_CHECK_EQUAL(res[2], 2);

  res = evaluate_coordinates_in_grid(res, grid);
  BOOST_CHECK_EQUAL(res[0], 2);
  BOOST_CHECK_EQUAL(res[1], 9);
  BOOST_CHECK_EQUAL(res[2], 8);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_critical_filtration_numerical_limits, T, list_of_tested_variants)
{
  BOOST_CHECK(std::numeric_limits<One_critical_filtration<T> >::has_infinity);
  BOOST_CHECK(std::numeric_limits<One_critical_filtration<T> >::infinity().is_inf());
  BOOST_CHECK(std::numeric_limits<One_critical_filtration<T> >::quiet_NaN().is_nan());
  BOOST_CHECK_THROW(std::numeric_limits<One_critical_filtration<T> >::max(), std::logic_error);
  auto max = std::numeric_limits<One_critical_filtration<T> >::max(3);
  BOOST_CHECK_EQUAL(max[0], std::numeric_limits<T>::max());
  BOOST_CHECK_EQUAL(max[1], std::numeric_limits<T>::max());
  BOOST_CHECK_EQUAL(max[2], std::numeric_limits<T>::max());
}

