/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cmath>        //std::isnan
#include <cstddef>      //std::size_t
#include <cstdint>      //std::int32_t
#include <limits>       //std::numerical_limits
#include <stdexcept>    //std::logic_error, std::out_of_range
#include <utility>      //std::swap, std::move
#include <vector>
#include <initializer_list>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_filtration"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Debug_utils.h>
#include <gudhi/Degree_rips_bifiltration.h>
#include <gudhi/Simplex_tree/filtration_value_utils.h>
#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>

using Gudhi::multi_filtration::Degree_rips_bifiltration;
using Gudhi::multi_filtration::Multi_parameter_filtration;
using Gudhi::multi_filtration::Dynamic_multi_parameter_filtration;

typedef boost::mpl::list<double, float, int> list_of_tested_variants;

template <class F, typename T, class F_alt>
void test_constructors(){
  F f0;
  BOOST_CHECK_EQUAL(f0.num_entries(), 2);
  BOOST_CHECK_EQUAL(f0.num_generators(), 1);
  BOOST_CHECK_EQUAL(f0.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f0(0,0), 0);
  BOOST_CHECK_EQUAL(f0(0,1), -F::T_inf);

  F f01(Gudhi::simplex_tree::empty_filtration_value_t{});
  BOOST_CHECK_EQUAL(f01.num_entries(), 0);
  BOOST_CHECK_EQUAL(f01.num_generators(), 0);
  BOOST_CHECK_EQUAL(f01.num_parameters(), 2);
  GUDHI_CHECK_code(BOOST_CHECK_THROW((f01[{0,0}]), std::out_of_range));

  F f1(3);
  BOOST_CHECK_EQUAL(f1.num_entries(), 2);
  BOOST_CHECK_EQUAL(f1.num_generators(), 1);
  BOOST_CHECK_EQUAL(f1.num_parameters(), 2);
  BOOST_CHECK_EQUAL((f1[{0,0}]), 0);
  BOOST_CHECK_EQUAL((f1[{0,1}]), -F::T_inf);
  
  F f2(3, 0);
  BOOST_CHECK_EQUAL(f2.num_entries(), 2);
  BOOST_CHECK_EQUAL(f2.num_generators(), 1);
  BOOST_CHECK_EQUAL(f2.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f2(0,0), 0);
  BOOST_CHECK_EQUAL(f2(0,1), 0);

  F f3({0,1,2});
  BOOST_CHECK_EQUAL(f3.num_entries(), 2);
  BOOST_CHECK_EQUAL(f3.num_generators(), 1);
  BOOST_CHECK_EQUAL(f3.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f3(0,0), 0);
  BOOST_CHECK_EQUAL(f3(0,1), 1);

  std::vector<T> v{0, 1, 1, 3, 2, 5};
  F f4(v.begin(), v.end());
  BOOST_CHECK_EQUAL(f4.num_entries(), 2);
  BOOST_CHECK_EQUAL(f4.num_generators(), 1);
  BOOST_CHECK_EQUAL(f4.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f4(0,0), 0);
  BOOST_CHECK_EQUAL(f4(0,1), 1);

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_THROW(F f9(v.begin(), v.end(), 3), std::logic_error);
    if constexpr (std::is_same_v<std::vector<T>, typename F::Underlying_container>){
      BOOST_CHECK_THROW(F f5(v, 3), std::logic_error);
      BOOST_CHECK_THROW(F f6(std::move(v), 3), std::logic_error);
    }
  } else {
    F f9(v.begin(), v.end(), 3);
    BOOST_CHECK_EQUAL(f9.num_entries(), 6);
    BOOST_CHECK_EQUAL(f9.num_generators(), 3);
    BOOST_CHECK_EQUAL(f9.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f9(0, 0), 0);
    BOOST_CHECK_EQUAL(f9(0, 1), 1);
    BOOST_CHECK_EQUAL(f9(1, 0), 1);
    BOOST_CHECK_EQUAL(f9(1, 1), 3);
    BOOST_CHECK_EQUAL(f9(2, 0), 2);
    BOOST_CHECK_EQUAL(f9(2, 1), 5);

    if constexpr (std::is_same_v<std::vector<T>, typename F::Underlying_container>){
      F f5(v, 3);
      BOOST_CHECK_EQUAL(f5.num_entries(), 12);
      BOOST_CHECK_EQUAL(f5.num_generators(), 6);
      BOOST_CHECK_EQUAL(f5.num_parameters(), 2);
      BOOST_CHECK_EQUAL((f5[{0,0}]), 0);
      BOOST_CHECK_EQUAL((f5[{0,1}]), 0);
      BOOST_CHECK_EQUAL((f5[{1,0}]), 1);
      BOOST_CHECK_EQUAL((f5[{1,1}]), 1);
      BOOST_CHECK_EQUAL((f5[{2,0}]), 2);
      BOOST_CHECK_EQUAL((f5[{2,1}]), 1);
      BOOST_CHECK_EQUAL((f5[{3,0}]), 3);
      BOOST_CHECK_EQUAL((f5[{3,1}]), 3);
      BOOST_CHECK_EQUAL((f5[{4,0}]), 4);
      BOOST_CHECK_EQUAL((f5[{4,1}]), 2);
      BOOST_CHECK_EQUAL((f5[{5,0}]), 5);
      BOOST_CHECK_EQUAL((f5[{5,1}]), 5);
    
      F f6(std::move(v), 3);
      BOOST_CHECK(v.empty());
      BOOST_CHECK_EQUAL(f6.num_entries(), 12);
      BOOST_CHECK_EQUAL(f6.num_generators(), 6);
      BOOST_CHECK_EQUAL(f6.num_parameters(), 2);
      BOOST_CHECK_EQUAL(f6(0,0), 0);
      BOOST_CHECK_EQUAL(f6(0,1), 0);
      BOOST_CHECK_EQUAL(f6(1,0), 1);
      BOOST_CHECK_EQUAL(f6(1,1), 1);
      BOOST_CHECK_EQUAL(f6(2,0), 2);
      BOOST_CHECK_EQUAL(f6(2,1), 1);
      BOOST_CHECK_EQUAL(f6(3,0), 3);
      BOOST_CHECK_EQUAL(f6(3,1), 3);
      BOOST_CHECK_EQUAL(f6(4,0), 4);
      BOOST_CHECK_EQUAL(f6(4,1), 2);
      BOOST_CHECK_EQUAL(f6(5,0), 5);
      BOOST_CHECK_EQUAL(f6(5,1), 5);
    }
  
    F f7(f9);
    BOOST_CHECK_EQUAL(f7.num_entries(), 6);
    BOOST_CHECK_EQUAL(f7.num_generators(), 3);
    BOOST_CHECK_EQUAL(f7.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f7(0, 0), 0);
    BOOST_CHECK_EQUAL(f7(0, 1), 1);
    BOOST_CHECK_EQUAL(f7(1, 0), 1);
    BOOST_CHECK_EQUAL(f7(1, 1), 3);
    BOOST_CHECK_EQUAL(f7(2, 0), 2);
    BOOST_CHECK_EQUAL(f7(2, 1), 5);
  
    F f8(std::move(f9));
    BOOST_CHECK_EQUAL(f8.num_entries(), 6);
    BOOST_CHECK_EQUAL(f8.num_generators(), 3);
    BOOST_CHECK_EQUAL(f8.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f8(0, 0), 0);
    BOOST_CHECK_EQUAL(f8(0, 1), 1);
    BOOST_CHECK_EQUAL(f8(1, 0), 1);
    BOOST_CHECK_EQUAL(f8(1, 1), 3);
    BOOST_CHECK_EQUAL(f8(2, 0), 2);
    BOOST_CHECK_EQUAL(f8(2, 1), 5);
  
    swap(f0, f8);
    BOOST_CHECK_EQUAL(f8.num_entries(), 2);
    BOOST_CHECK_EQUAL(f8.num_generators(), 1);
    BOOST_CHECK_EQUAL(f8.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f8(0,0), 0);
    BOOST_CHECK_EQUAL(f8(0,1), -F::T_inf);
    BOOST_CHECK_EQUAL(f0.num_entries(), 6);
    BOOST_CHECK_EQUAL(f0.num_generators(), 3);
    BOOST_CHECK_EQUAL(f0.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f0(0,0), 0);
    BOOST_CHECK_EQUAL(f0(0,1), 1);
    BOOST_CHECK_EQUAL(f0(1,0), 1);
    BOOST_CHECK_EQUAL(f0(1,1), 3);
    BOOST_CHECK_EQUAL(f0(2,0), 2);
    BOOST_CHECK_EQUAL(f0(2,1), 5);
  }

  F_alt f10({0,1,2});
  F f11(f10);
  BOOST_CHECK_EQUAL(f11.num_entries(), 2);
  BOOST_CHECK_EQUAL(f11.num_generators(), 1);
  BOOST_CHECK_EQUAL(f11.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f11(0,0), 0);
  BOOST_CHECK_EQUAL(f11(0,1), 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_constructors, T, list_of_tested_variants)
{
  test_constructors<Degree_rips_bifiltration<T>, T, Degree_rips_bifiltration<std::int32_t> >();
  test_constructors<Degree_rips_bifiltration<T, false, true>, T, Degree_rips_bifiltration<std::int32_t> >();
}

template <class F, typename T, class F_alt>
void test_utilities(){
  F f0({0, 1, 2});
  bool test = std::is_same_v<decltype(f0(0, 0)), T&>;
  BOOST_CHECK(test);

  F_alt f2 = f0.template as_type<float>();
  test = std::is_same_v<decltype(f2(0, 0)), float&>;
  BOOST_CHECK(test);
  BOOST_CHECK_EQUAL(f2.num_generators(), 1);
  BOOST_CHECK_EQUAL(f2.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f2(0, 0), 0.);
  BOOST_CHECK_EQUAL(f2(0, 1), 1.);

  BOOST_CHECK(!f0.is_plus_inf());
  BOOST_CHECK(!f0.is_minus_inf());
  BOOST_CHECK(!f0.is_nan());
  BOOST_CHECK(f0.is_finite());

  F f3;
  BOOST_CHECK(!f3.is_plus_inf());
  BOOST_CHECK(f3.is_minus_inf());
  BOOST_CHECK(!f3.is_nan());
  BOOST_CHECK(!f3.is_finite());

  F f4 = F::minus_inf(3);
  BOOST_CHECK(!f3.is_plus_inf());
  BOOST_CHECK(f3.is_minus_inf());
  BOOST_CHECK(!f3.is_nan());
  BOOST_CHECK(!f3.is_finite());

  F f5 = F::inf(3);
  BOOST_CHECK(f5.is_plus_inf());
  BOOST_CHECK(!f5.is_minus_inf());
  BOOST_CHECK(!f5.is_nan());
  BOOST_CHECK(!f5.is_finite());

  if constexpr (std::numeric_limits<F>::has_quiet_NaN){
    F f6 = F::nan(3);
    BOOST_CHECK(!f6.is_plus_inf());
    BOOST_CHECK(!f6.is_minus_inf());
    BOOST_CHECK(f6.is_nan());
    BOOST_CHECK(!f6.is_finite());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_utilities, T, list_of_tested_variants)
{
  test_utilities<Degree_rips_bifiltration<T>, T, Degree_rips_bifiltration<float> >();
  test_utilities<Degree_rips_bifiltration<T, false, true>, T, Degree_rips_bifiltration<float> >();
}

template <class F, typename T>
void test_comparators(){
  const int num_param = 3;
  std::vector<T> v1, v2, v3, v4;

  if constexpr (F::ensures_1_criticality()) {
    v1 = {0, 1};
    v2 = {0, 0};
    v3 = {0, 2};
    v4 = {0, 1};
  } else {
    v1 = {0, 4, 1, 3, 2, 2};
    v2 = {0, 0, 1, -2, 2, -1, 3, 0};
    v3 = {0, 5};
    v4 = {0, 5, 1, 3, 2, 1};
  }

  F f1(v1.begin(), v1.end(), num_param);
  F f2(v2.begin(), v2.end(), num_param);
  F f3(v3.begin(), v3.end(), num_param);
  F f4(v4.begin(), v4.end(), num_param);

  BOOST_CHECK(!(f1 < f1));
  BOOST_CHECK(!(f1 < f2));
  BOOST_CHECK(f1 < f3);
  BOOST_CHECK(!(f1 < f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 < F::nan(num_param)));
  BOOST_CHECK(f1 < F::inf(num_param));
  BOOST_CHECK(!(f1 < F::minus_inf(num_param)));

  BOOST_CHECK(f1 <= f1);
  BOOST_CHECK(!(f1 <= f2));
  BOOST_CHECK(f1 <= f3);
  if constexpr (F::ensures_1_criticality()) BOOST_CHECK(f1 <= f4);
  else BOOST_CHECK(!(f1 <= f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 <= F::nan(num_param)));
  BOOST_CHECK(f1 <= F::inf(num_param));
  BOOST_CHECK(!(f1 <= F::minus_inf(num_param)));

  BOOST_CHECK(!(f1 > f1));
  BOOST_CHECK(f1 > f2);
  BOOST_CHECK(!(f1 > f3));
  BOOST_CHECK(!(f1 > f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 > F::nan(num_param)));
  BOOST_CHECK(!(f1 > F::inf(num_param)));
  BOOST_CHECK(f1 > F::minus_inf(num_param));

  BOOST_CHECK(f1 >= f1);
  BOOST_CHECK(f1 >= f2);
  BOOST_CHECK(!(f1 >= f3));
  if constexpr (F::ensures_1_criticality()) BOOST_CHECK(f1 >= f4);
  else BOOST_CHECK(!(f1 >= f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 >= F::nan(num_param)));
  BOOST_CHECK(!(f1 >= F::inf(num_param)));
  BOOST_CHECK(f1 >= F::minus_inf(num_param));

  BOOST_CHECK(f1 == f1);
  BOOST_CHECK(!(f1 == f2));
  BOOST_CHECK(!(f1 == f3));
  if constexpr (F::ensures_1_criticality()) BOOST_CHECK(f1 == f4);
  else BOOST_CHECK(!(f1 == f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 == F::nan(num_param)));
  BOOST_CHECK(!(f1 == F::inf(num_param)));
  BOOST_CHECK(!(f1 == F::minus_inf(num_param)));

  BOOST_CHECK(!(f1 != f1));
  BOOST_CHECK(f1 != f2);
  BOOST_CHECK(f1 != f3);
  if constexpr (F::ensures_1_criticality()) BOOST_CHECK(!(f1 != f4));
  else BOOST_CHECK(f1 != f4);
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(f1 != F::nan(num_param));
  BOOST_CHECK(f1 != F::inf(num_param));
  BOOST_CHECK(f1 != F::minus_inf(num_param));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_comparators, T, list_of_tested_variants)
{
  test_comparators<Degree_rips_bifiltration<T>, T>();
  test_comparators<Degree_rips_bifiltration<T, false, true>, T>();
}

template <class F, typename T>
void test_operators(){
  const int num_param = 3;

  F f({0, -1});
  F f2({0, 2});
  F f3({0, -F::T_inf});
  F f4({0, F::T_inf});
  // TODO: tests with more than 1 generator

  F res = -f;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 1);
  BOOST_CHECK((-F::inf(num_param)).is_minus_inf());
  BOOST_CHECK((-F::minus_inf(num_param)).is_plus_inf());
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK((-F::nan(num_param)).is_nan());

  res = f - f2;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), -3);

  res = f - f3;
  BOOST_CHECK_EQUAL(res(0,0), f4(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f4(0,1));

  res = f3 - f;
  BOOST_CHECK_EQUAL(res(0,0), f3(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f3(0,1));

  res = 5 - f;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 6);

  res = f - 5;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), -6);
  BOOST_CHECK((f - F::inf(num_param)).is_minus_inf());
  BOOST_CHECK((F::inf(num_param) - f).is_plus_inf());
  BOOST_CHECK((f - F::minus_inf(num_param)).is_plus_inf());
  BOOST_CHECK((F::minus_inf(num_param) - f).is_minus_inf());
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK((f - F::nan(num_param)) == f);
    BOOST_CHECK((F::nan(num_param) - f).is_nan());
  }

  res = f3 - f3;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(std::isnan(res(0,1)));
  else BOOST_CHECK_EQUAL(res(0,1), 0);
  res = f3 - f4;
  BOOST_CHECK_EQUAL(res(0,0), f3(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f3(0,1));

  res = f + f2;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 1);

  res = f + f3;
  BOOST_CHECK_EQUAL(res(0,0), f3(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f3(0,1));

  res = 5 + f;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 4);

  res = f + 5;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 4);

  BOOST_CHECK((f + F::inf(num_param)).is_plus_inf());
  BOOST_CHECK((F::inf(num_param) + f).is_plus_inf());
  BOOST_CHECK((f + F::minus_inf(num_param)).is_minus_inf());
  BOOST_CHECK((F::minus_inf(num_param) + f).is_minus_inf());
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK((f + F::nan(num_param)) == f);
    BOOST_CHECK((F::nan(num_param) + f).is_nan());
  }

  res = f3 + f4;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(std::isnan(res(0,1)));
  else BOOST_CHECK_EQUAL(res(0,1), 0);
  res = f3 + f3;
  BOOST_CHECK_EQUAL(res(0,0), f3(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f3(0,1));

  res = f * f2;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), -2);

  res = f * f3;
  BOOST_CHECK_EQUAL(res(0,0), f4(0,0));
  BOOST_CHECK_EQUAL(res(0,0), -f4(0,0));

  res = 0 * f3;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }

  res = 5 * f;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), -5);

  res = f * 5;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), -5);

  res = f * F::inf(num_param);
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), -F::T_inf);

  res = F::inf(num_param) * f;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), -F::T_inf);

  res = f * F::minus_inf(num_param);
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), F::T_inf);

  res = F::minus_inf(num_param) * f;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), F::T_inf);

  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK((f * F::nan(num_param)) == f);
    BOOST_CHECK((F::nan(num_param) * f).is_nan());
  }

  res = f3 * f3;
  BOOST_CHECK(res.is_plus_inf());
  res = f3 * f4;
  BOOST_CHECK(res.is_minus_inf());

  res = f / f2;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), static_cast<T>(-0.5));

  res = f / f3;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 0);

  res = f3 / f;
  BOOST_CHECK_EQUAL(res(0,0), f4(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f4(0,1));

  res = f3 / 0;
  BOOST_CHECK_EQUAL(res(0,0), f4(0,0));
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }

  res = 5 / f;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), -5);

  res = f / 5;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), static_cast<T>(-0.2));

  res = f / F::inf(num_param);
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 0);
  res = F::inf(num_param) / f;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), f3(0,1));

  res = f / F::minus_inf(num_param);
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 0);
  res = F::minus_inf(num_param) / f;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), f4(0,1));

  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK((f / F::nan(num_param)) == f);
    BOOST_CHECK((F::nan(num_param) / f).is_nan());
  }

  res = f3 / f3;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(std::isnan(res(0,1)));
  else BOOST_CHECK_EQUAL(res(0,1), 0);
  res = f3 / f4;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(std::isnan(res(0,1)));
  else BOOST_CHECK_EQUAL(res(0,1), 0);
  res = f / F({0, 0, 0});
  BOOST_CHECK_EQUAL(res(0,0), 0);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) BOOST_CHECK(std::isnan(res(0,1)));
  else BOOST_CHECK_EQUAL(res(0,1), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_operators, T, list_of_tested_variants)
{
  test_operators<Degree_rips_bifiltration<T>, T>();
  test_operators<Degree_rips_bifiltration<T, false, true>, T>();
}

template <class F, typename T>
void test_modifiers(){
  const int num_param = 3;
  std::vector<T> v;

  if constexpr (F::ensures_1_criticality()) {
    v = {0, 1};
  } else {
    v = {0, 1, 1, 0, 2, 2};
  }

  F f(v.begin(), v.end(), num_param);
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 1);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 1);
    BOOST_CHECK_EQUAL(f(1, 0), 1);
    BOOST_CHECK_EQUAL(f(1, 1), 0);
    BOOST_CHECK_EQUAL(f(2, 0), 2);
    BOOST_CHECK_EQUAL(f(2, 1), 2);
  }

  f.push_to_least_common_upper_bound({0, 1});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 1);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 1);
    BOOST_CHECK_EQUAL(f(1, 0), 1);
    BOOST_CHECK_EQUAL(f(1, 1), 1);
    BOOST_CHECK_EQUAL(f(2, 0), 2);
    BOOST_CHECK_EQUAL(f(2, 1), 2);
  }

  f.push_to_least_common_upper_bound({0,3});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 3);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 3);
    BOOST_CHECK_EQUAL(f(1, 0), 1);
    BOOST_CHECK_EQUAL(f(1, 1), 3);
    BOOST_CHECK_EQUAL(f(2, 0), 2);
    BOOST_CHECK_EQUAL(f(2, 1), 3);
  }

  f.push_to_least_common_upper_bound(F::minus_inf(num_param));
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 3);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 3);
    BOOST_CHECK_EQUAL(f(1, 0), 1);
    BOOST_CHECK_EQUAL(f(1, 1), 3);
    BOOST_CHECK_EQUAL(f(2, 0), 2);
    BOOST_CHECK_EQUAL(f(2, 1), 3);
  }

  f.push_to_least_common_upper_bound(F::inf(num_param));
  BOOST_CHECK(f.is_plus_inf());

  if constexpr (std::numeric_limits<F>::has_quiet_NaN){
    f.push_to_least_common_upper_bound(F::nan(num_param));
    BOOST_CHECK(f.is_plus_inf());
  }

  f.pull_to_greatest_common_lower_bound({5, 4});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 4);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 4);
    BOOST_CHECK_EQUAL(f(1, 0), 1);
    BOOST_CHECK_EQUAL(f(1, 1), 4);
    BOOST_CHECK_EQUAL(f(2, 0), 2);
    BOOST_CHECK_EQUAL(f(2, 1), 4);
  }

  f.pull_to_greatest_common_lower_bound({5,-1});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), -1);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), -1);
    BOOST_CHECK_EQUAL(f(1, 0), 1);
    BOOST_CHECK_EQUAL(f(1, 1), -1);
    BOOST_CHECK_EQUAL(f(2, 0), 2);
    BOOST_CHECK_EQUAL(f(2, 1), -1);
  }

  f.pull_to_greatest_common_lower_bound(F::inf(num_param));
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), -1);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), -1);
    BOOST_CHECK_EQUAL(f(1, 0), 1);
    BOOST_CHECK_EQUAL(f(1, 1), -1);
    BOOST_CHECK_EQUAL(f(2, 0), 2);
    BOOST_CHECK_EQUAL(f(2, 1), -1);
  }

  f.pull_to_greatest_common_lower_bound(F::minus_inf(num_param));
  BOOST_CHECK(f.is_minus_inf());

  if constexpr (std::numeric_limits<F>::has_quiet_NaN){
    f.pull_to_greatest_common_lower_bound(F::nan(num_param));
    BOOST_CHECK(f.is_minus_inf());
  }

  std::vector<std::vector<int> > grid = {{0, 1, 2, 3}, {0, 3, 6, 9}, {0, 4, 8, 16}};

  f.push_to_least_common_upper_bound({0, 7});
  f.project_onto_grid(grid, true);
  BOOST_CHECK_EQUAL(f(0, 0), 0);
  BOOST_CHECK_EQUAL(f(0, 1), 3);

  f.push_to_least_common_upper_bound({0, 7});
  f.project_onto_grid(grid, false);
  BOOST_CHECK_EQUAL(f(0, 0), 0);
  BOOST_CHECK_EQUAL(f(0, 1), 9);

  if constexpr (!F::ensures_1_criticality()) {
    f.set_num_generators(5);
    BOOST_CHECK_EQUAL(f.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f.num_generators(), 5);
    BOOST_CHECK_EQUAL(f.num_entries(), 10);
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 9);
    BOOST_CHECK_EQUAL(f(1, 0), 1);
    BOOST_CHECK_EQUAL(f(1, 1), 9);
    BOOST_CHECK_EQUAL(f(2, 0), 2);
    BOOST_CHECK_EQUAL(f(2, 1), 9);
    BOOST_CHECK_EQUAL(f(3, 0), 3);
    BOOST_CHECK_EQUAL(f(3, 1), -F::T_inf);
    BOOST_CHECK_EQUAL(f(4, 0), 4);
    BOOST_CHECK_EQUAL(f(4, 1), -F::T_inf);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_modifiers, T, list_of_tested_variants)
{
  test_modifiers<Degree_rips_bifiltration<T>, T>();
  test_modifiers<Degree_rips_bifiltration<T, false, true>, T>();
}

template <class F, typename T>
void test_add_generators(){
  const int num_param = 2;

  F f({0, 1});
  BOOST_CHECK_EQUAL(f.num_generators(), 1);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), 0);
  BOOST_CHECK_EQUAL(f(0, 1), 1);

  bool res = f.add_generator({1, 2});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), 0);
  BOOST_CHECK_EQUAL(f(0, 1), 1);
  BOOST_CHECK_EQUAL(f(1, 0), 1);
  BOOST_CHECK_EQUAL(f(1, 1), 2);

  res = f.add_generator({0, -2});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), 0);
  BOOST_CHECK_EQUAL(f(0, 1), -2);
  BOOST_CHECK_EQUAL(f(1, 0), 1);
  BOOST_CHECK_EQUAL(f(1, 1), 2);

  res = f.add_generator({0, 3});
  BOOST_CHECK(!res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), 0);
  BOOST_CHECK_EQUAL(f(0, 1), -2);
  BOOST_CHECK_EQUAL(f(1, 0), 1);
  BOOST_CHECK_EQUAL(f(1, 1), 2);

  res = f.add_generator({0, F::T_inf});
  BOOST_CHECK(!res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), 0);
  BOOST_CHECK_EQUAL(f(0, 1), -2);
  BOOST_CHECK_EQUAL(f(1, 0), 1);
  BOOST_CHECK_EQUAL(f(1, 1), 2);

  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    res = f.add_generator({0, std::numeric_limits<T>::quiet_NaN()});
    BOOST_CHECK(!res);
    BOOST_CHECK_EQUAL(f.num_generators(), 2);
    BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), -2);
    BOOST_CHECK_EQUAL(f(1, 0), 1);
    BOOST_CHECK_EQUAL(f(1, 1), 2);
  }

  res = f.add_generator({0, -F::T_inf});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), 0);
  BOOST_CHECK_EQUAL(f(0, 1), -F::T_inf);
  BOOST_CHECK_EQUAL(f(1, 0), 1);
  BOOST_CHECK_EQUAL(f(1, 1), 2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_add_generators, T, list_of_tested_variants)
{
  test_add_generators<Degree_rips_bifiltration<T>, T>();

  Degree_rips_bifiltration<T, false, true> f({0, 1});
  BOOST_CHECK_THROW(f.add_generator({1, 1}), std::logic_error);
}

template <class F, typename T>
void test_friends(){
  F f({1, 2}, 2);

  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(6))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, F({4, 5, 3}, 2)), static_cast<T>(std::sqrt(T(2))));
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {2, 3, 5, 9}), 3);
  F ff = factorize_below(f);
  BOOST_CHECK(ff == F({0, 1}));
  BOOST_CHECK(ff <= f);
  ff = factorize_above(f);
  BOOST_CHECK(ff == F({0, 2}));
  BOOST_CHECK(ff >= f);

  f.add_guaranteed_generator({2, 0});
  BOOST_CHECK_EQUAL(f.num_generators(), 3);
  BOOST_CHECK_EQUAL(f.num_parameters(), 2);

  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(10))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, F({4, 5, 3}, 2)),
                    static_cast<T>(std::sqrt(T(2))));
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {2, 3, 5, 9}), 3);
  ff = factorize_below(f);
  BOOST_CHECK(ff == F({0, 0}));
  BOOST_CHECK(ff <= f);
  ff = factorize_above(f);
  BOOST_CHECK(ff == F({0, 2}));
  BOOST_CHECK(ff >= f);

  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    T nan = std::numeric_limits<T>::quiet_NaN();
    std::vector<T> v = {0, nan, 1, 2, 2, nan};
    F f2(v.begin(), v.end(), 3);

    BOOST_CHECK(std::isnan(compute_norm(f2)));
    BOOST_CHECK(std::isnan(compute_euclidean_distance_to(f2, std::initializer_list<T>{0, 2})));
    BOOST_CHECK(std::isnan(compute_linear_projection(f2, {0, 3})));
    F f2f = factorize_below(f2);
    BOOST_CHECK_EQUAL(f2f(0, 0), 0);
    BOOST_CHECK_EQUAL(f2f(0, 1), 2);
    f2f = factorize_above(f2);
    BOOST_CHECK_EQUAL(f2f(0, 0), 0);
    BOOST_CHECK_EQUAL(f2f(0, 1), 2);
  }

  f(0,1) = 1;
  f(1,1) = 7;
  f(2,1) = 5;

  std::vector<std::vector<int> > grid = {{0, 1, 2, 3}, {0, 3, 6, 9}, {0, 4, 8, 16}};
  auto res = compute_coordinates_in_grid(f, grid);
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res.num_generators(), 3);
  BOOST_CHECK_EQUAL(f.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f.num_generators(), 3);
  BOOST_CHECK_EQUAL(res(0, 1), 1);
  BOOST_CHECK_EQUAL(res(1, 1), 3);
  BOOST_CHECK_EQUAL(res(2, 1), 2);

  res = evaluate_coordinates_in_grid(res, grid);
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res.num_generators(), 3);
  BOOST_CHECK_EQUAL(res(0, 1), 3);
  BOOST_CHECK_EQUAL(res(1, 1), 9);
  BOOST_CHECK_EQUAL(res(2, 1), 6);
}

template <class F, typename T>
void test_friends_1_critical(){
  F f({0, 1});

  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(1))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, F({0,3})), 2);
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {2, 3, 5, 9}), 3);
  BOOST_CHECK(factorize_below(f) == f);
  BOOST_CHECK(factorize_above(f) == f);

  f(0,1) = 7;

  std::vector<std::vector<int> > grid = {{0, 1, 2, 3}, {0, 3, 6, 9}, {0, 4, 8, 16}};
  auto res = compute_coordinates_in_grid(f, grid);
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res.num_generators(), 1);
    BOOST_CHECK_EQUAL(f.num_generators(), 1);
  BOOST_CHECK_EQUAL(res(0, 0), 0);
  BOOST_CHECK_EQUAL(res(0, 1), 3);

  res = evaluate_coordinates_in_grid(res, grid);
  BOOST_CHECK_EQUAL(res.num_generators(), 1);
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res(0, 0), 0);
  BOOST_CHECK_EQUAL(res(0, 1), 9);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_friends, T, list_of_tested_variants)
{
  test_friends<Degree_rips_bifiltration<T>, T>();
  test_friends_1_critical<Degree_rips_bifiltration<T, false, true>, T>();
}

template <class F, typename T>
void test_unify_intersect(){
  const int num_param = 2;

  std::vector<T> v1 = {0,4,1,3,2,1,3,2};
  F f1(v1.begin(), v1.end(), num_param);

  std::vector<T> v2 = {0,5,1,2,2,1};
  F f2(v2.begin(), v2.end(), num_param);

  bool modified = unify_lifetimes(f1, f2);
  BOOST_CHECK(modified);
  BOOST_CHECK(f1.num_parameters() == num_param);
  BOOST_CHECK(f1.num_generators() == 4);
  BOOST_CHECK_EQUAL(f1(0,0), 0);
  BOOST_CHECK_EQUAL(f1(0,1), 4);
  BOOST_CHECK_EQUAL(f1(1,0), 1);
  BOOST_CHECK_EQUAL(f1(1,1), 2);
  BOOST_CHECK_EQUAL(f1(2,0), 2);
  BOOST_CHECK_EQUAL(f1(2,1), 1);
  BOOST_CHECK_EQUAL(f1(3,0), 3);
  BOOST_CHECK_EQUAL(f1(3,1), 2);

  std::vector<T> v3 = {0,4,1,3,2,1,3,2};
  F f3(v3.begin(), v3.end(), num_param);

  modified = intersect_lifetimes(f3, f2);
  BOOST_CHECK(modified);
  BOOST_CHECK(f3.num_parameters() == num_param);
  BOOST_CHECK(f3.num_generators() == 4);
  BOOST_CHECK_EQUAL(f3(0,0), 0);
  BOOST_CHECK_EQUAL(f3(0,1), 5);
  BOOST_CHECK_EQUAL(f3(1,0), 1);
  BOOST_CHECK_EQUAL(f3(1,1), 3);
  BOOST_CHECK_EQUAL(f3(2,0), 2);
  BOOST_CHECK_EQUAL(f3(2,1), 1);
  BOOST_CHECK_EQUAL(f3(3,0), 3);
  BOOST_CHECK_EQUAL(f3(3,1), 1);
}

template <class F, typename T>
void test_unify_intersect_1_critical(){
  const int num_param = 2;

  std::vector<T> v1 = {0,5};
  F f1(v1.begin(), v1.end(), num_param);

  std::vector<T> v2 = {0,8};
  F f2(v2.begin(), v2.end(), num_param);

  bool modified = unify_lifetimes(f1, f2);
  BOOST_CHECK(!modified);
  BOOST_CHECK(f1.num_parameters() == num_param);
  BOOST_CHECK(f1.num_generators() == 1);
  BOOST_CHECK_EQUAL(f1(0,0), 0);
  BOOST_CHECK_EQUAL(f1(0,1), 5);

  modified = intersect_lifetimes(f1, f2);
  BOOST_CHECK(modified);
  BOOST_CHECK(f1.num_parameters() == num_param);
  BOOST_CHECK(f1.num_generators() == 1);
  BOOST_CHECK_EQUAL(f1(0,0), 0);
  BOOST_CHECK_EQUAL(f1(0,1), 8);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_unify_intersect, T, list_of_tested_variants)
{
  test_unify_intersect<Degree_rips_bifiltration<T>, T>();
  test_unify_intersect_1_critical<Degree_rips_bifiltration<T, false, true>, T>();
}

template <class F, typename T>
void test_serialize(){
  std::vector<T> v = {0,5,1,3,2,2};
  F f(v.begin(), v.end(), 2);
  BOOST_CHECK(f.num_parameters() == 2);
  BOOST_CHECK(f.num_generators() == 3);
  BOOST_CHECK_EQUAL(f(0,0), 0);
  BOOST_CHECK_EQUAL(f(0,1), 5);
  BOOST_CHECK_EQUAL(f(1,0), 1);
  BOOST_CHECK_EQUAL(f(1,1), 3);
  BOOST_CHECK_EQUAL(f(2,0), 2);
  BOOST_CHECK_EQUAL(f(2,1), 2);

  char* buffer = new char[256];
  std::size_t serializationSize = get_serialization_size_of(f);
  
  char* ptr = buffer;
  ptr = serialize_value_to_char_buffer(f, ptr);
  BOOST_CHECK_EQUAL(static_cast<std::size_t>(ptr - buffer), serializationSize);

  const char* c_ptr = buffer;
  F f3(Gudhi::simplex_tree::empty_filtration_value_t{});
  c_ptr = deserialize_value_from_char_buffer(f3, c_ptr);
  BOOST_CHECK_EQUAL(static_cast<std::size_t>(c_ptr - buffer), serializationSize);
  BOOST_CHECK(f3.num_parameters() == 2);
  BOOST_CHECK(f3.num_generators() == 3);
  BOOST_CHECK_EQUAL(f3(0,0), 0);
  BOOST_CHECK_EQUAL(f3(0,1), 5);
  BOOST_CHECK_EQUAL(f3(1,0), 1);
  BOOST_CHECK_EQUAL(f3(1,1), 3);
  BOOST_CHECK_EQUAL(f3(2,0), 2);
  BOOST_CHECK_EQUAL(f3(2,1), 2);

  delete[] buffer;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_serialize, T, list_of_tested_variants)
{
  test_serialize<Degree_rips_bifiltration<T>, T>();
}

template <class F, typename T>
void test_co(){
  F f;
  BOOST_CHECK(f.num_parameters() == 2);
  BOOST_CHECK(f.num_generators() == 1);
  BOOST_CHECK_EQUAL(f(0,0), 0);
  BOOST_CHECK_EQUAL(f(0,1), F::T_inf);

  BOOST_CHECK_THROW(F::inf(), std::logic_error);
  BOOST_CHECK(!f.is_plus_inf());
  BOOST_CHECK(!f.is_minus_inf());
  BOOST_CHECK(!f.is_nan());
  BOOST_CHECK(f.is_finite());

  F f6 = F::minus_inf(3);
  bool change = f6.add_generator({0, F::T_inf});
  BOOST_CHECK(change);
  BOOST_CHECK_EQUAL(f6(0,0), 0);
  BOOST_CHECK_EQUAL(f6(0,1), F::T_inf);

  if constexpr (F::ensures_1_criticality()) {
    std::vector<T> v = {0, 1};
    F f2(v.begin(), v.end(), 3);
    BOOST_CHECK_EQUAL(compute_linear_projection(f2, {2,3,5,9}), 3);
  } else {
    std::vector<T> v = {0, 1, 1, 2, 2, 4};
    F f2(v.begin(), v.end(), 3);
    BOOST_CHECK_EQUAL(compute_linear_projection(f2, {2,3,5,9}), 16);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_co, T, list_of_tested_variants)
{
  test_co<Degree_rips_bifiltration<T, true>, T>();
  test_co<Degree_rips_bifiltration<T, true, true>, T>();
}

template <class F, typename T, bool Co>
void test_numerical_limits(){
  const int num_param = 3;

  if constexpr (Co) BOOST_CHECK(!std::numeric_limits<F>::has_infinity);
  else BOOST_CHECK(std::numeric_limits<F>::has_infinity);
  BOOST_CHECK(std::numeric_limits<F>::has_quiet_NaN);
  
  BOOST_CHECK(std::numeric_limits<F>::quiet_NaN().is_nan());
  BOOST_CHECK(std::numeric_limits<F>::minus_infinity().is_minus_inf());
  if constexpr (Co) {
    BOOST_CHECK_THROW(std::numeric_limits<F>::infinity(), std::logic_error);
    BOOST_CHECK_THROW(std::numeric_limits<F>::max(), std::logic_error);
  } else {
    BOOST_CHECK(std::numeric_limits<F>::infinity().is_plus_inf());
    auto max = std::numeric_limits<F>::max();
    BOOST_CHECK_EQUAL(max(0,0), 0);
    BOOST_CHECK_EQUAL(max(0,1), std::numeric_limits<T>::max());
  }

  BOOST_CHECK(std::numeric_limits<F>::quiet_NaN(num_param).is_nan());
  BOOST_CHECK(std::numeric_limits<F>::minus_infinity(num_param).is_minus_inf());
  if constexpr (Co) {
    BOOST_CHECK_THROW(std::numeric_limits<F>::infinity(num_param), std::logic_error);
    BOOST_CHECK_THROW(std::numeric_limits<F>::max(num_param), std::logic_error);
  } else {
    BOOST_CHECK(std::numeric_limits<F>::infinity().is_plus_inf());
    auto max = std::numeric_limits<F>::max(num_param);
    BOOST_CHECK_EQUAL(max(0,0), 0);
    BOOST_CHECK_EQUAL(max(0,1), std::numeric_limits<T>::max());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_numerical_limits, T, list_of_tested_variants)
{
  test_numerical_limits<Degree_rips_bifiltration<T>, T, false>();
  test_numerical_limits<Degree_rips_bifiltration<T, false, true>, T, false>();
  test_numerical_limits<Degree_rips_bifiltration<T, true>, T, true>();
  test_numerical_limits<Degree_rips_bifiltration<T, true, true>, T, true>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_converters, T, list_of_tested_variants)
{
  std::vector<T> v = {5, 6, 3, 4};
  Degree_rips_bifiltration<T> f(std::move(v), 2);
  BOOST_CHECK(f.num_parameters() == 2);
  BOOST_CHECK(f.num_generators() == 4);
  BOOST_CHECK_EQUAL(f(0,0), 0);
  BOOST_CHECK_EQUAL(f(0,1), 5);
  BOOST_CHECK_EQUAL(f(1,0), 1);
  BOOST_CHECK_EQUAL(f(1,1), 6);
  BOOST_CHECK_EQUAL(f(2,0), 2);
  BOOST_CHECK_EQUAL(f(2,1), 3);
  BOOST_CHECK_EQUAL(f(3,0), 3);
  BOOST_CHECK_EQUAL(f(3,1), 4);

  Multi_parameter_filtration<T> f11 = f.convert_to_multi_parameter_filtration();
  BOOST_CHECK(f11.num_parameters() == 2);
  BOOST_CHECK(f11.num_generators() == 2);
  BOOST_CHECK_EQUAL(f11(0,0), 0);
  BOOST_CHECK_EQUAL(f11(0,1), 5);
  BOOST_CHECK_EQUAL(f11(1,0), 2);
  BOOST_CHECK_EQUAL(f11(1,1), 3);

  Multi_parameter_filtration<T> f12 = f.convert_to_non_simplified_multi_parameter_filtration();
  BOOST_CHECK(f12.num_parameters() == 2);
  BOOST_CHECK(f12.num_generators() == 4);
  BOOST_CHECK_EQUAL(f12(0,0), 0);
  BOOST_CHECK_EQUAL(f12(0,1), 5);
  BOOST_CHECK_EQUAL(f12(1,0), 1);
  BOOST_CHECK_EQUAL(f12(1,1), 6);
  BOOST_CHECK_EQUAL(f12(2,0), 2);
  BOOST_CHECK_EQUAL(f12(2,1), 3);
  BOOST_CHECK_EQUAL(f12(3,0), 3);
  BOOST_CHECK_EQUAL(f12(3,1), 4);

  Dynamic_multi_parameter_filtration<T> f21 = f.convert_to_dynamic_multi_parameter_filtration();
  BOOST_CHECK(f21.num_parameters() == 2);
  BOOST_CHECK(f21.num_generators() == 2);
  BOOST_CHECK_EQUAL(f21(0,0), 0);
  BOOST_CHECK_EQUAL(f21(0,1), 5);
  BOOST_CHECK_EQUAL(f21(1,0), 2);
  BOOST_CHECK_EQUAL(f21(1,1), 3);

  Dynamic_multi_parameter_filtration<T> f22 = f.convert_to_non_simplified_dynamic_multi_parameter_filtration();
  BOOST_CHECK(f22.num_parameters() == 2);
  BOOST_CHECK(f22.num_generators() == 4);
  BOOST_CHECK_EQUAL(f22(0,0), 0);
  BOOST_CHECK_EQUAL(f22(0,1), 5);
  BOOST_CHECK_EQUAL(f22(1,0), 1);
  BOOST_CHECK_EQUAL(f22(1,1), 6);
  BOOST_CHECK_EQUAL(f22(2,0), 2);
  BOOST_CHECK_EQUAL(f22(2,1), 3);
  BOOST_CHECK_EQUAL(f22(3,0), 3);
  BOOST_CHECK_EQUAL(f22(3,1), 4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_converters_co, T, list_of_tested_variants)
{
  std::vector<T> v = {5, 6, 3, 4};
  Degree_rips_bifiltration<T, true> f(std::move(v), 2);
  BOOST_CHECK(f.num_parameters() == 2);
  BOOST_CHECK(f.num_generators() == 4);
  BOOST_CHECK_EQUAL(f(0,0), 0);
  BOOST_CHECK_EQUAL(f(0,1), 5);
  BOOST_CHECK_EQUAL(f(1,0), 1);
  BOOST_CHECK_EQUAL(f(1,1), 6);
  BOOST_CHECK_EQUAL(f(2,0), 2);
  BOOST_CHECK_EQUAL(f(2,1), 3);
  BOOST_CHECK_EQUAL(f(3,0), 3);
  BOOST_CHECK_EQUAL(f(3,1), 4);

  Multi_parameter_filtration<T, true> f11 = f.convert_to_multi_parameter_filtration();
  BOOST_CHECK(f11.num_parameters() == 2);
  BOOST_CHECK(f11.num_generators() == 2);
  BOOST_CHECK_EQUAL(f11(0,0), 1);
  BOOST_CHECK_EQUAL(f11(0,1), 6);
  BOOST_CHECK_EQUAL(f11(1,0), 3);
  BOOST_CHECK_EQUAL(f11(1,1), 4);

  Multi_parameter_filtration<T, true> f12 = f.convert_to_non_simplified_multi_parameter_filtration();
  BOOST_CHECK(f12.num_parameters() == 2);
  BOOST_CHECK(f12.num_generators() == 4);
  BOOST_CHECK_EQUAL(f12(0,0), 0);
  BOOST_CHECK_EQUAL(f12(0,1), 5);
  BOOST_CHECK_EQUAL(f12(1,0), 1);
  BOOST_CHECK_EQUAL(f12(1,1), 6);
  BOOST_CHECK_EQUAL(f12(2,0), 2);
  BOOST_CHECK_EQUAL(f12(2,1), 3);
  BOOST_CHECK_EQUAL(f12(3,0), 3);
  BOOST_CHECK_EQUAL(f12(3,1), 4);

  Dynamic_multi_parameter_filtration<T, true> f21 = f.convert_to_dynamic_multi_parameter_filtration();
  BOOST_CHECK(f21.num_parameters() == 2);
  BOOST_CHECK(f21.num_generators() == 2);
  BOOST_CHECK_EQUAL(f21(0,0), 1);
  BOOST_CHECK_EQUAL(f21(0,1), 6);
  BOOST_CHECK_EQUAL(f21(1,0), 3);
  BOOST_CHECK_EQUAL(f21(1,1), 4);

  Dynamic_multi_parameter_filtration<T, true> f22 = f.convert_to_non_simplified_dynamic_multi_parameter_filtration();
  BOOST_CHECK(f22.num_parameters() == 2);
  BOOST_CHECK(f22.num_generators() == 4);
  BOOST_CHECK_EQUAL(f22(0,0), 0);
  BOOST_CHECK_EQUAL(f22(0,1), 5);
  BOOST_CHECK_EQUAL(f22(1,0), 1);
  BOOST_CHECK_EQUAL(f22(1,1), 6);
  BOOST_CHECK_EQUAL(f22(2,0), 2);
  BOOST_CHECK_EQUAL(f22(2,1), 3);
  BOOST_CHECK_EQUAL(f22(3,0), 3);
  BOOST_CHECK_EQUAL(f22(3,1), 4);
}

