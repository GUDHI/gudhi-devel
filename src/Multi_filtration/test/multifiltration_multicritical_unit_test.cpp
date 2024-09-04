/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cmath>
#include <limits>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_filtration"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_critical_filtration.h>

using Gudhi::multi_filtration::Multi_critical_filtration;

typedef boost::mpl::list<double, float, int> list_of_tested_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_constructors, T, list_of_tested_variants)
{
  Multi_critical_filtration<T, false> f01;
  BOOST_CHECK(f01.num_parameters() == 1);
  BOOST_CHECK(f01.num_generators() == 1);
  BOOST_CHECK(f01[0][0] == -Multi_critical_filtration<T>::Generator::T_inf);
  Multi_critical_filtration<T, true> f02;
  BOOST_CHECK(f02.num_parameters() == 1);
  BOOST_CHECK(f02.num_generators() == 1);
  BOOST_CHECK(f02[0][0] == Multi_critical_filtration<T>::Generator::T_inf);

  Multi_critical_filtration<T> f1(3);
  BOOST_CHECK(f1.num_parameters() == 3);
  BOOST_CHECK(f1.num_generators() == 1);
  BOOST_CHECK(f1[0][0] == -Multi_critical_filtration<T>::Generator::T_inf);
  BOOST_CHECK(f1[0][1] == -Multi_critical_filtration<T>::Generator::T_inf);
  BOOST_CHECK(f1[0][2] == -Multi_critical_filtration<T>::Generator::T_inf);

  Multi_critical_filtration<T> f2(3, 0);
  BOOST_CHECK(f2.num_parameters() == 3);
  BOOST_CHECK(f2.num_generators() == 1);
  BOOST_CHECK(f2[0][0] == 0);
  BOOST_CHECK(f2[0][1] == 0);
  BOOST_CHECK(f2[0][2] == 0);

  Multi_critical_filtration<T> f3({0, 1, 2});
  BOOST_CHECK(f3.num_parameters() == 3);
  BOOST_CHECK(f3.num_generators() == 1);
  BOOST_CHECK(f3[0][0] == 0);
  BOOST_CHECK(f3[0][1] == 1);
  BOOST_CHECK(f3[0][2] == 2);

  std::vector<T> v{0, 1, 2};
  Multi_critical_filtration<T> f41(v);
  BOOST_CHECK(f41.num_parameters() == 3);
  BOOST_CHECK(f41.num_generators() == 1);
  BOOST_CHECK(f41[0][0] == 0);
  BOOST_CHECK(f41[0][1] == 1);
  BOOST_CHECK(f41[0][2] == 2);

  Multi_critical_filtration<T> f5(v.begin(), v.end());
  BOOST_CHECK(f5.num_generators() == 1);
  BOOST_CHECK(f5.num_parameters() == 3);
  BOOST_CHECK(f5[0][0] == 0);
  BOOST_CHECK(f5[0][1] == 1);
  BOOST_CHECK(f5[0][2] == 2);

  Multi_critical_filtration<T> f42(std::move(v));
  BOOST_CHECK(f42.num_parameters() == 3);
  BOOST_CHECK(f42.num_generators() == 1);
  BOOST_CHECK(f42[0][0] == 0);
  BOOST_CHECK(f42[0][1] == 1);
  BOOST_CHECK(f42[0][2] == 2);

  std::vector<typename Multi_critical_filtration<T>::Generator> v2{{0, 1, 2}, {3, 4, 5}};
  Multi_critical_filtration<T> f6(std::move(v2));
  BOOST_CHECK(f6.num_parameters() == 3);
  BOOST_CHECK(f6.num_generators() == 2);
  BOOST_CHECK(f6[0][0] == 0);
  BOOST_CHECK(f6[0][1] == 1);
  BOOST_CHECK(f6[0][2] == 2);
  BOOST_CHECK(f6[1][0] == 3);
  BOOST_CHECK(f6[1][1] == 4);
  BOOST_CHECK(f6[1][2] == 5);

  Multi_critical_filtration<T> f7(f6);
  BOOST_CHECK(f7.num_parameters() == 3);
  BOOST_CHECK(f7.num_generators() == 2);
  BOOST_CHECK(f7[0][0] == 0);
  BOOST_CHECK(f7[0][1] == 1);
  BOOST_CHECK(f7[0][2] == 2);
  BOOST_CHECK(f7[1][0] == 3);
  BOOST_CHECK(f7[1][1] == 4);
  BOOST_CHECK(f7[1][2] == 5);

  f01 = f5;
  BOOST_CHECK(f01.num_generators() == 1);
  BOOST_CHECK(f01.num_parameters() == 3);
  BOOST_CHECK(f01[0][0] == 0);
  BOOST_CHECK(f01[0][1] == 1);
  BOOST_CHECK(f01[0][2] == 2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_utilities, T, list_of_tested_variants)
{
  Multi_critical_filtration<T> f({0, 1, 2});
  bool test = std::is_same_v<decltype(f[0][0]), T&>;
  BOOST_CHECK(test);

  Multi_critical_filtration<float> f2 = f.template as_type<float>();
  test = std::is_same_v<decltype(f2[0][0]), float&>;
  BOOST_CHECK(test);
  BOOST_CHECK(f2.num_generators() == 1);
  BOOST_CHECK(f2.num_parameters() == 3);
  BOOST_CHECK(f2[0][0] == 0.);
  BOOST_CHECK(f2[0][1] == 1.);
  BOOST_CHECK(f2[0][2] == 2.);

  BOOST_CHECK(!f.is_inf());
  BOOST_CHECK(!f.is_minus_inf());
  BOOST_CHECK(!f.is_nan());
  BOOST_CHECK(f.is_finite());

  Multi_critical_filtration<T, false> f31;
  BOOST_CHECK(!f31.is_inf());
  BOOST_CHECK(f31.is_minus_inf());
  BOOST_CHECK(!f31.is_nan());
  BOOST_CHECK(!f31.is_finite());

  Multi_critical_filtration<T, true> f32;
  BOOST_CHECK(f32.is_inf());
  BOOST_CHECK(!f32.is_minus_inf());
  BOOST_CHECK(!f32.is_nan());
  BOOST_CHECK(!f32.is_finite());

  //{-inf, -inf, -inf} is considered finite as the user is supposed to updates the values to something else
  //the idea is just to reserve space and to give the possibility to use `f4[i] =`
  //if the value should really be -inf, use `f4(1)` or `f4 = minus_inf()` instead.
  Multi_critical_filtration<T> f4(3);
  BOOST_CHECK(!f4.is_inf());
  BOOST_CHECK(!f4.is_minus_inf());
  BOOST_CHECK(!f4.is_nan());
  BOOST_CHECK(f4.is_finite());

  Multi_critical_filtration<T> f5(1);
  BOOST_CHECK(!f5.is_inf());
  BOOST_CHECK(f5.is_minus_inf());
  BOOST_CHECK(!f5.is_nan());
  BOOST_CHECK(!f5.is_finite());

  auto f6 = Multi_critical_filtration<T>::nan();
  BOOST_CHECK(!f6.is_inf());
  BOOST_CHECK(!f6.is_minus_inf());
  BOOST_CHECK(f6.is_nan());
  BOOST_CHECK(!f6.is_finite());

  auto f7 = Multi_critical_filtration<T>::inf();
  BOOST_CHECK(f7.is_inf());
  BOOST_CHECK(!f7.is_minus_inf());
  BOOST_CHECK(!f7.is_nan());
  BOOST_CHECK(!f7.is_finite());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_comparators, T, list_of_tested_variants)
{
  Multi_critical_filtration<T> f1({{0, 1, 2}, {-1, 4, 5}});
  Multi_critical_filtration<T> f2({{-2, 0, 1}, {-5, 0, 0}, {-2, -5, -1}});
  Multi_critical_filtration<T> f3({4, 5, 6});
  Multi_critical_filtration<T> f4({{0, 0, 1}, {-4, 5, 6}});

  BOOST_CHECK(!(f1 < f1));
  BOOST_CHECK(!(f1 < f2));
  BOOST_CHECK(f1 < f3);
  BOOST_CHECK(!(f1 < f4));
  BOOST_CHECK(!(f1 < Multi_critical_filtration<T>::nan()));
  BOOST_CHECK(f1 < Multi_critical_filtration<T>::inf());
  BOOST_CHECK(!(f1 < Multi_critical_filtration<T>::minus_inf()));

  BOOST_CHECK(f1 <= f1);
  BOOST_CHECK(!(f1 <= f2));
  BOOST_CHECK(f1 <= f3);
  BOOST_CHECK(!(f1 <= f4));
  BOOST_CHECK(!(f1 <= Multi_critical_filtration<T>::nan()));
  BOOST_CHECK(f1 <= Multi_critical_filtration<T>::inf());
  BOOST_CHECK(!(f1 <= Multi_critical_filtration<T>::minus_inf()));

  BOOST_CHECK(!(f1 > f1));
  BOOST_CHECK(f1 > f2);
  BOOST_CHECK(!(f1 > f3));
  BOOST_CHECK(!(f1 > f4));
  BOOST_CHECK(!(f1 > Multi_critical_filtration<T>::nan()));
  BOOST_CHECK(!(f1 > Multi_critical_filtration<T>::inf()));
  BOOST_CHECK(f1 > Multi_critical_filtration<T>::minus_inf());

  BOOST_CHECK(f1 >= f1);
  BOOST_CHECK(f1 >= f2);
  BOOST_CHECK(!(f1 >= f3));
  BOOST_CHECK(!(f1 >= f4));
  BOOST_CHECK(!(f1 >= Multi_critical_filtration<T>::nan()));
  BOOST_CHECK(!(f1 >= Multi_critical_filtration<T>::inf()));
  BOOST_CHECK(f1 >= Multi_critical_filtration<T>::minus_inf());

  BOOST_CHECK(f1 == f1);
  BOOST_CHECK(!(f1 == f2));
  BOOST_CHECK(!(f1 == f3));
  BOOST_CHECK(!(f1 == f4));
  BOOST_CHECK(!(f1 == Multi_critical_filtration<T>::nan()));
  BOOST_CHECK(!(f1 == Multi_critical_filtration<T>::inf()));
  BOOST_CHECK(!(f1 == Multi_critical_filtration<T>::minus_inf()));

  BOOST_CHECK(!(f1 != f1));
  BOOST_CHECK(f1 != f2);
  BOOST_CHECK(f1 != f3);
  BOOST_CHECK(f1 != f4);
  BOOST_CHECK(f1 != Multi_critical_filtration<T>::nan());
  BOOST_CHECK(f1 != Multi_critical_filtration<T>::inf());
  BOOST_CHECK(f1 != Multi_critical_filtration<T>::minus_inf());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_modifiers, T, list_of_tested_variants)
{
  Multi_critical_filtration<T> f({{0, 1, 2}, {-3, 1, 7}});
  BOOST_CHECK_EQUAL(f[0][0], 0);
  BOOST_CHECK_EQUAL(f[0][1], 1);
  BOOST_CHECK_EQUAL(f[0][2], 2);
  BOOST_CHECK_EQUAL(f[1][0], -3);
  BOOST_CHECK_EQUAL(f[1][1], 1);
  BOOST_CHECK_EQUAL(f[1][2], 7);

  f.push_to_least_common_upper_bound({-1, 5, 6});
  BOOST_CHECK_EQUAL(f[0][0], 0);
  BOOST_CHECK_EQUAL(f[0][1], 5);
  BOOST_CHECK_EQUAL(f[0][2], 6);
  BOOST_CHECK_EQUAL(f[1][0], -1);
  BOOST_CHECK_EQUAL(f[1][1], 5);
  BOOST_CHECK_EQUAL(f[1][2], 7);

  f.push_to_least_common_upper_bound({-1, -5, -6});
  BOOST_CHECK_EQUAL(f[0][0], 0);
  BOOST_CHECK_EQUAL(f[0][1], 5);
  BOOST_CHECK_EQUAL(f[0][2], 6);
  BOOST_CHECK_EQUAL(f[1][0], -1);
  BOOST_CHECK_EQUAL(f[1][1], 5);
  BOOST_CHECK_EQUAL(f[1][2], 7);

  f.push_to_least_common_upper_bound(Multi_critical_filtration<T>::Generator::minus_inf());
  BOOST_CHECK_EQUAL(f[0][0], 0);
  BOOST_CHECK_EQUAL(f[0][1], 5);
  BOOST_CHECK_EQUAL(f[0][2], 6);
  BOOST_CHECK_EQUAL(f[1][0], -1);
  BOOST_CHECK_EQUAL(f[1][1], 5);
  BOOST_CHECK_EQUAL(f[1][2], 7);

  f.push_to_least_common_upper_bound(Multi_critical_filtration<T>::Generator::inf());
  BOOST_CHECK(f.is_inf());

  f.push_to_least_common_upper_bound(Multi_critical_filtration<T>::Generator::nan());
  BOOST_CHECK(f.is_inf());

  f.pull_to_greatest_common_lower_bound({-1, 5, 6});
  BOOST_CHECK_EQUAL(f[0][0], -1);
  BOOST_CHECK_EQUAL(f[0][1], 5);
  BOOST_CHECK_EQUAL(f[0][2], 6);

  f.pull_to_greatest_common_lower_bound({1, 8, 9});
  BOOST_CHECK_EQUAL(f[0][0], -1);
  BOOST_CHECK_EQUAL(f[0][1], 5);
  BOOST_CHECK_EQUAL(f[0][2], 6);

  f.pull_to_greatest_common_lower_bound(Multi_critical_filtration<T>::Generator::inf());
  BOOST_CHECK_EQUAL(f[0][0], -1);
  BOOST_CHECK_EQUAL(f[0][1], 5);
  BOOST_CHECK_EQUAL(f[0][2], 6);

  f.pull_to_greatest_common_lower_bound(Multi_critical_filtration<T>::Generator::minus_inf());
  BOOST_CHECK(f.is_minus_inf());

  f.pull_to_greatest_common_lower_bound(Multi_critical_filtration<T>::Generator::nan());
  BOOST_CHECK(f.is_minus_inf());

  std::vector<std::vector<int> > grid = {{0, 2, 4, 8}, {0, 3, 6, 9}, {0, 4, 8, 16}};

  f.push_to_least_common_upper_bound({1, 7, 5});
  f.project_onto_grid(grid, true);
  BOOST_CHECK_EQUAL(f[0][0], 1);
  BOOST_CHECK_EQUAL(f[0][1], 3);
  BOOST_CHECK_EQUAL(f[0][2], 2);

  f.push_to_least_common_upper_bound({1, 7, 5});
  f.project_onto_grid(grid, false);
  BOOST_CHECK_EQUAL(f[0][0], 2);
  BOOST_CHECK_EQUAL(f[0][1], 9);
  BOOST_CHECK_EQUAL(f[0][2], 8);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_add, T, list_of_tested_variants)
{
  Multi_critical_filtration<T> f({0, 1, 2});
  BOOST_CHECK_EQUAL(f.num_generators(), 1);
  BOOST_CHECK_EQUAL(f.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f[0][0], 0);
  BOOST_CHECK_EQUAL(f[0][1], 1);
  BOOST_CHECK_EQUAL(f[0][2], 2);

  bool res = f.add_generator({-3, 1, 7});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f[0][0], 0);
  BOOST_CHECK_EQUAL(f[0][1], 1);
  BOOST_CHECK_EQUAL(f[0][2], 2);
  BOOST_CHECK_EQUAL(f[1][0], -3);
  BOOST_CHECK_EQUAL(f[1][1], 1);
  BOOST_CHECK_EQUAL(f[1][2], 7);

  res = f.add_generator({-1, -2, -3});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f[0][0], -3);
  BOOST_CHECK_EQUAL(f[0][1], 1);
  BOOST_CHECK_EQUAL(f[0][2], 7);
  BOOST_CHECK_EQUAL(f[1][0], -1);
  BOOST_CHECK_EQUAL(f[1][1], -2);
  BOOST_CHECK_EQUAL(f[1][2], -3);

  res = f.add_generator({8, 9, 10});
  BOOST_CHECK(!res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f[0][0], -3);
  BOOST_CHECK_EQUAL(f[0][1], 1);
  BOOST_CHECK_EQUAL(f[0][2], 7);
  BOOST_CHECK_EQUAL(f[1][0], -1);
  BOOST_CHECK_EQUAL(f[1][1], -2);
  BOOST_CHECK_EQUAL(f[1][2], -3);

  res = f.add_generator(Multi_critical_filtration<T>::Generator::inf());
  BOOST_CHECK(!res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f[0][0], -3);
  BOOST_CHECK_EQUAL(f[0][1], 1);
  BOOST_CHECK_EQUAL(f[0][2], 7);
  BOOST_CHECK_EQUAL(f[1][0], -1);
  BOOST_CHECK_EQUAL(f[1][1], -2);
  BOOST_CHECK_EQUAL(f[1][2], -3);

  res = f.add_generator(Multi_critical_filtration<T>::Generator::nan());
  BOOST_CHECK(!res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f[0][0], -3);
  BOOST_CHECK_EQUAL(f[0][1], 1);
  BOOST_CHECK_EQUAL(f[0][2], 7);
  BOOST_CHECK_EQUAL(f[1][0], -1);
  BOOST_CHECK_EQUAL(f[1][1], -2);
  BOOST_CHECK_EQUAL(f[1][2], -3);

  res = f.add_generator(Multi_critical_filtration<T>::Generator::minus_inf());
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 1);
  BOOST_CHECK_EQUAL(f.num_parameters(), 1);
  BOOST_CHECK_EQUAL(f[0][0], -Multi_critical_filtration<T>::Generator::T_inf);

  std::vector<typename Multi_critical_filtration<T>::Generator> v{
      {0, 1, 2},
      typename Multi_critical_filtration<T>::Generator(0),
      Multi_critical_filtration<T>::Generator::inf(),
      {0, 1, 2},
      Multi_critical_filtration<T>::Generator::nan(),
      typename Multi_critical_filtration<T>::Generator(0),
      Multi_critical_filtration<T>::Generator::minus_inf()};

  Multi_critical_filtration<T> f2(v);
  f2.remove_empty_generators(false);
  BOOST_CHECK_EQUAL(f2[0][0], 0);
  BOOST_CHECK_EQUAL(f2[0][1], 1);
  BOOST_CHECK_EQUAL(f2[0][2], 2);
  BOOST_CHECK(f2[1].is_inf());
  BOOST_CHECK_EQUAL(f2[2][0], 0);
  BOOST_CHECK_EQUAL(f2[2][1], 1);
  BOOST_CHECK_EQUAL(f2[2][2], 2);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(f2.num_generators(), 5);
    BOOST_CHECK(f2[3].is_nan());
    BOOST_CHECK(f2[4].is_minus_inf());
  } else {
    BOOST_CHECK_EQUAL(f2.num_generators(), 4);
    BOOST_CHECK(f2[3].is_minus_inf());
  }

  Multi_critical_filtration<T> f3(v);
  f3.remove_empty_generators(true);
  BOOST_CHECK_EQUAL(f3[0][0], 0);
  BOOST_CHECK_EQUAL(f3[0][1], 1);
  BOOST_CHECK_EQUAL(f3[0][2], 2);
  BOOST_CHECK_EQUAL(f3[1][0], 0);
  BOOST_CHECK_EQUAL(f3[1][1], 1);
  BOOST_CHECK_EQUAL(f3[1][2], 2);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(f3.num_generators(), 3);
    BOOST_CHECK(f3[2].is_nan());
  } else {
    BOOST_CHECK_EQUAL(f3.num_generators(), 2);
  }

  Multi_critical_filtration<T> f4(v);
  f4.simplify();
  BOOST_CHECK_EQUAL(f4.num_generators(), 1);
  BOOST_CHECK(f4[0].is_minus_inf());

  v.pop_back();
  Multi_critical_filtration<T> f5(v);
  f5.simplify();
  BOOST_CHECK_EQUAL(f5.num_generators(), 1);
  BOOST_CHECK_EQUAL(f5[0][0], 0);
  BOOST_CHECK_EQUAL(f5[0][1], 1);
  BOOST_CHECK_EQUAL(f5[0][2], 2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_friends, T, list_of_tested_variants)
{
  Multi_critical_filtration<T> f({{0, 1, 2}, {2, 0, 4}});
  
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {2,3,5,9}), 13);
  BOOST_CHECK(factorize_below(f) == typename Multi_critical_filtration<T>::Generator({0, 0, 2}));
  BOOST_CHECK(factorize_above(f) == typename Multi_critical_filtration<T>::Generator({2, 1, 4}));

  Multi_critical_filtration<T, true> f2({{0, 1, 2}, {2, 0, 4}});
  BOOST_CHECK_EQUAL(compute_linear_projection(f2, {2,3,5,9}), 24);

  f[0] = {1, 7, 5};

  std::vector<std::vector<int> > grid = {{0, 2, 4, 8}, {0, 3, 6, 9}, {0, 4, 8, 16}};
  auto res = compute_coordinates_in_grid(f, grid);
  BOOST_CHECK_EQUAL(res[0][0], 1);
  BOOST_CHECK_EQUAL(res[0][1], 3);
  BOOST_CHECK_EQUAL(res[0][2], 2);
  BOOST_CHECK_EQUAL(res[1][0], 1);
  BOOST_CHECK_EQUAL(res[1][1], 0);
  BOOST_CHECK_EQUAL(res[1][2], 1);

  res = evaluate_coordinates_in_grid(res, grid);
  BOOST_CHECK_EQUAL(res[0][0], 2);
  BOOST_CHECK_EQUAL(res[0][1], 0);
  BOOST_CHECK_EQUAL(res[0][2], 4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_numerical_limits, T, list_of_tested_variants)
{
  BOOST_CHECK(std::numeric_limits<Multi_critical_filtration<T> >::has_infinity);
  BOOST_CHECK(std::numeric_limits<Multi_critical_filtration<T> >::infinity().is_inf());
  BOOST_CHECK(std::numeric_limits<Multi_critical_filtration<T> >::quiet_NaN().is_nan());
  BOOST_CHECK_THROW(std::numeric_limits<Multi_critical_filtration<T> >::max(), std::logic_error);
  auto max = std::numeric_limits<Multi_critical_filtration<T> >::max(2, 3);
  BOOST_CHECK_EQUAL(max[0][0], std::numeric_limits<T>::max());
  BOOST_CHECK_EQUAL(max[0][1], std::numeric_limits<T>::max());
  BOOST_CHECK_EQUAL(max[0][2], std::numeric_limits<T>::max());
  BOOST_CHECK_EQUAL(max[1][0], std::numeric_limits<T>::max());
  BOOST_CHECK_EQUAL(max[1][1], std::numeric_limits<T>::max());
  BOOST_CHECK_EQUAL(max[1][2], std::numeric_limits<T>::max());
}

