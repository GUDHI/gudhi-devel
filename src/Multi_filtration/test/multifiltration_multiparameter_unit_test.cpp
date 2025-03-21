/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cmath>    //std::isnan
#include <cstdint>  //std::int32_t
#include <limits>   //std::numerical_limits
#include <type_traits>
#include <utility>  //std::swap, std::move
#include <vector>
#include <initializer_list>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_filtration"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>

using Gudhi::multi_filtration::Multi_parameter_filtration;
using Gudhi::multi_filtration::Dynamic_multi_parameter_filtration;

typedef boost::mpl::list<double, float, int> list_of_tested_variants;

template <class F, typename T, class F_alt>
void test_constructors(){
  F f0;
  BOOST_CHECK_EQUAL(f0.num_entries(), 2);
  BOOST_CHECK_EQUAL(f0.num_generators(), 1);
  BOOST_CHECK_EQUAL(f0.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f0(0,0), -F::T_inf);
  BOOST_CHECK_EQUAL(f0(0,1), -F::T_inf);
  
  F f1(3);
  BOOST_CHECK_EQUAL(f1.num_entries(), 3);
  BOOST_CHECK_EQUAL(f1.num_generators(), 1);
  BOOST_CHECK_EQUAL(f1.num_parameters(), 3);
  BOOST_CHECK_EQUAL((f1[{0,0}]), -F::T_inf);
  BOOST_CHECK_EQUAL((f1[{0,1}]), -F::T_inf);
  BOOST_CHECK_EQUAL((f1[{0,2}]), -F::T_inf);

  F f2(3, 0);
  BOOST_CHECK_EQUAL(f2.num_entries(), 3);
  BOOST_CHECK_EQUAL(f2.num_generators(), 1);
  BOOST_CHECK_EQUAL(f2.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f2(0,0), 0);
  BOOST_CHECK_EQUAL(f2(0,1), 0);
  BOOST_CHECK_EQUAL(f2(0,2), 0);

  F f3({0,1,2});
  BOOST_CHECK_EQUAL(f3.num_entries(), 3);
  BOOST_CHECK_EQUAL(f3.num_generators(), 1);
  BOOST_CHECK_EQUAL(f3.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f3(0,0), 0);
  BOOST_CHECK_EQUAL(f3(0,1), 1);
  BOOST_CHECK_EQUAL(f3(0,2), 2);

  std::vector<T> v{0, 1, 2, 3, 4, 5};
  F f4(v.begin(), v.end());
  BOOST_CHECK_EQUAL(f4.num_entries(), 6);
  BOOST_CHECK_EQUAL(f4.num_generators(), 1);
  BOOST_CHECK_EQUAL(f4.num_parameters(), 6);
  BOOST_CHECK_EQUAL(f4(0,0), 0);
  BOOST_CHECK_EQUAL(f4(0,1), 1);
  BOOST_CHECK_EQUAL(f4(0,2), 2);
  BOOST_CHECK_EQUAL(f4(0,3), 3);
  BOOST_CHECK_EQUAL(f4(0,4), 4);
  BOOST_CHECK_EQUAL(f4(0,5), 5);

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_THROW(F f9(v.begin(), v.end(), 3), std::logic_error);
    if constexpr (std::is_same_v<std::vector<T>, typename F::Underlying_container>){
      BOOST_CHECK_THROW(F f5(v, 3), std::logic_error);
      BOOST_CHECK_THROW(F f6(std::move(v), 3), std::logic_error);
    }
  } else {
    F f9(v.begin(), v.end(), 3);
    BOOST_CHECK_EQUAL(f9.num_entries(), 6);
    BOOST_CHECK_EQUAL(f9.num_generators(), 2);
    BOOST_CHECK_EQUAL(f9.num_parameters(), 3);
    BOOST_CHECK_EQUAL(f9(0, 0), 0);
    BOOST_CHECK_EQUAL(f9(0, 1), 1);
    BOOST_CHECK_EQUAL(f9(0, 2), 2);
    BOOST_CHECK_EQUAL(f9(1, 0), 3);
    BOOST_CHECK_EQUAL(f9(1, 1), 4);
    BOOST_CHECK_EQUAL(f9(1, 2), 5);

    if constexpr (std::is_same_v<std::vector<T>, typename F::Underlying_container>){
      F f5(v, 3);
      BOOST_CHECK_EQUAL(f5.num_entries(), 6);
      BOOST_CHECK_EQUAL(f5.num_generators(), 2);
      BOOST_CHECK_EQUAL(f5.num_parameters(), 3);
      BOOST_CHECK_EQUAL((f5[{0,0}]), 0);
      BOOST_CHECK_EQUAL((f5[{0,1}]), 1);
      BOOST_CHECK_EQUAL((f5[{0,2}]), 2);
      BOOST_CHECK_EQUAL((f5[{1,0}]), 3);
      BOOST_CHECK_EQUAL((f5[{1,1}]), 4);
      BOOST_CHECK_EQUAL((f5[{1,2}]), 5);
    
      F f6(std::move(v), 3);
      BOOST_CHECK(v.empty());
      BOOST_CHECK_EQUAL(f6.num_entries(), 6);
      BOOST_CHECK_EQUAL(f6.num_generators(), 2);
      BOOST_CHECK_EQUAL(f6.num_parameters(), 3);
      BOOST_CHECK_EQUAL(f6(0,0), 0);
      BOOST_CHECK_EQUAL(f6(0,1), 1);
      BOOST_CHECK_EQUAL(f6(0,2), 2);
      BOOST_CHECK_EQUAL(f6(1,0), 3);
      BOOST_CHECK_EQUAL(f6(1,1), 4);
      BOOST_CHECK_EQUAL(f6(1,2), 5);
    }
  
    F f7(f9);
    BOOST_CHECK_EQUAL(f7.num_entries(), 6);
    BOOST_CHECK_EQUAL(f7.num_generators(), 2);
    BOOST_CHECK_EQUAL(f7.num_parameters(), 3);
    BOOST_CHECK_EQUAL(f7(0,0), 0);
    BOOST_CHECK_EQUAL(f7(0,1), 1);
    BOOST_CHECK_EQUAL(f7(0,2), 2);
    BOOST_CHECK_EQUAL(f7(1,0), 3);
    BOOST_CHECK_EQUAL(f7(1,1), 4);
    BOOST_CHECK_EQUAL(f7(1,2), 5);
  
    F f8(std::move(f9));
    BOOST_CHECK_EQUAL(f8.num_entries(), 6);
    BOOST_CHECK_EQUAL(f8.num_generators(), 2);
    BOOST_CHECK_EQUAL(f8.num_parameters(), 3);
    BOOST_CHECK_EQUAL(f8(0,0), 0);
    BOOST_CHECK_EQUAL(f8(0,1), 1);
    BOOST_CHECK_EQUAL(f8(0,2), 2);
    BOOST_CHECK_EQUAL(f8(1,0), 3);
    BOOST_CHECK_EQUAL(f8(1,1), 4);
    BOOST_CHECK_EQUAL(f8(1,2), 5);
  
    swap(f0, f8);
    BOOST_CHECK_EQUAL(f8.num_entries(), 2);
    BOOST_CHECK_EQUAL(f8.num_generators(), 1);
    BOOST_CHECK_EQUAL(f8.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f8(0,0), -F::T_inf);
    BOOST_CHECK_EQUAL(f8(0,1), -F::T_inf);
    BOOST_CHECK_EQUAL(f0.num_entries(), 6);
    BOOST_CHECK_EQUAL(f0.num_generators(), 2);
    BOOST_CHECK_EQUAL(f0.num_parameters(), 3);
    BOOST_CHECK_EQUAL(f0(0,0), 0);
    BOOST_CHECK_EQUAL(f0(0,1), 1);
    BOOST_CHECK_EQUAL(f0(0,2), 2);
    BOOST_CHECK_EQUAL(f0(1,0), 3);
    BOOST_CHECK_EQUAL(f0(1,1), 4);
    BOOST_CHECK_EQUAL(f0(1,2), 5);
  }

  F_alt f10({0,1,2});
  F f11(f10);
  BOOST_CHECK_EQUAL(f11.num_entries(), 3);
  BOOST_CHECK_EQUAL(f11.num_generators(), 1);
  BOOST_CHECK_EQUAL(f11.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f11(0,0), 0);
  BOOST_CHECK_EQUAL(f11(0,1), 1);
  BOOST_CHECK_EQUAL(f11(0,2), 2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_constructors, T, list_of_tested_variants)
{
  test_constructors<Multi_parameter_filtration<T>, T, Multi_parameter_filtration<std::int32_t> >();
  test_constructors<Multi_parameter_filtration<T, false, true>, T, Multi_parameter_filtration<std::int32_t> >();

  test_constructors<Dynamic_multi_parameter_filtration<T>, T, Dynamic_multi_parameter_filtration<std::int32_t> >();
  test_constructors<Dynamic_multi_parameter_filtration<T, false, true>,
                    T,
                    Dynamic_multi_parameter_filtration<std::int32_t> >();
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
  BOOST_CHECK_EQUAL(f2.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f2(0, 0), 0.);
  BOOST_CHECK_EQUAL(f2(0, 1), 1.);
  BOOST_CHECK_EQUAL(f2(0, 2), 2.);

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

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_utilities, T, list_of_tested_variants)
{
  test_utilities<Multi_parameter_filtration<T>, T, Multi_parameter_filtration<float> >();
  test_utilities<Multi_parameter_filtration<T, false, true>, T, Multi_parameter_filtration<float> >();

  test_utilities<Dynamic_multi_parameter_filtration<T>, T, Dynamic_multi_parameter_filtration<float> >();
  test_utilities<Dynamic_multi_parameter_filtration<T, false, true>, T, Dynamic_multi_parameter_filtration<float> >();
}

template <class F, typename T>
void test_comparators(){
  const int num_param = 3;
  std::vector<T> v1, v2, v3, v4;

  if constexpr (F::ensures_1_criticality()) {
    v1 = {-1, 1, 2};
    v2 = {-2, 0, 1};
    v3 = {0, 2, 3};
    v4 = {5, -1, 2};
  } else {
    v1 = {-1, 4, 5, 0, 1, 2};
    v2 = {-5, 0, 0, -2, -5, -1, -2, 0, 1};
    v3 = {4, 5, 6};
    v4 = {-4, 5, 6, 0, 0, 1};
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
  BOOST_CHECK(!(f1 <= f4));
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
  BOOST_CHECK(!(f1 >= f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 >= F::nan(num_param)));
  BOOST_CHECK(!(f1 >= F::inf(num_param)));
  BOOST_CHECK(f1 >= F::minus_inf(num_param));

  BOOST_CHECK(f1 == f1);
  BOOST_CHECK(!(f1 == f2));
  BOOST_CHECK(!(f1 == f3));
  BOOST_CHECK(!(f1 == f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 == F::nan(num_param)));
  BOOST_CHECK(!(f1 == F::inf(num_param)));
  BOOST_CHECK(!(f1 == F::minus_inf(num_param)));

  BOOST_CHECK(!(f1 != f1));
  BOOST_CHECK(f1 != f2);
  BOOST_CHECK(f1 != f3);
  BOOST_CHECK(f1 != f4);
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(f1 != F::nan(num_param));
  BOOST_CHECK(f1 != F::inf(num_param));
  BOOST_CHECK(f1 != F::minus_inf(num_param));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_comparators, T, list_of_tested_variants)
{
  test_comparators<Multi_parameter_filtration<T>, T>();
  test_comparators<Multi_parameter_filtration<T, false, true>, T>();

  test_comparators<Dynamic_multi_parameter_filtration<T>, T>();
  test_comparators<Dynamic_multi_parameter_filtration<T, false, true>, T>();
}

template <class F, typename T>
void test_operators(){
  const int num_param = 3;

  F f({-10, 0, 1});
  F f2({5, 2, -1});
  F f3({-F::T_inf, F::T_inf, -F::T_inf});
  F f4({F::T_inf, -F::T_inf, F::T_inf});
  // TODO: tests with more than 1 generator

  F res = -f;
  BOOST_CHECK_EQUAL(res(0,0), 10);
  BOOST_CHECK_EQUAL(res(0,1), 0);
  BOOST_CHECK_EQUAL(res(0,2), -1);
  BOOST_CHECK((-F::inf(num_param)).is_minus_inf());
  BOOST_CHECK((-F::minus_inf(num_param)).is_plus_inf());
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK((-F::nan(num_param)).is_nan());

  res = f - f2;
  BOOST_CHECK_EQUAL(res(0,0), -15);
  BOOST_CHECK_EQUAL(res(0,1), -2);
  BOOST_CHECK_EQUAL(res(0,2), 2);

  res = f - f3;
  BOOST_CHECK_EQUAL(res(0,0), f4(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f4(0,1));
  BOOST_CHECK_EQUAL(res(0,2), f4(0,2));

  res = f3 - f;
  BOOST_CHECK_EQUAL(res(0,0), f3(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f3(0,1));
  BOOST_CHECK_EQUAL(res(0,2), f3(0,2));

  res = T(5) - f;
  BOOST_CHECK_EQUAL(res(0,0), 15);
  BOOST_CHECK_EQUAL(res(0,1), 5);
  BOOST_CHECK_EQUAL(res(0,2), 4);

  res = f - T(5);
  BOOST_CHECK_EQUAL(res(0,0), -15);
  BOOST_CHECK_EQUAL(res(0,1), -5);
  BOOST_CHECK_EQUAL(res(0,2), -4);
  BOOST_CHECK((f - F::inf(num_param)).is_minus_inf());
  BOOST_CHECK((F::inf(num_param) - f).is_plus_inf());
  BOOST_CHECK((f - F::minus_inf(num_param)).is_plus_inf());
  BOOST_CHECK((F::minus_inf(num_param) - f).is_minus_inf());
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK((f - F::nan(num_param)).is_nan());
    BOOST_CHECK((F::nan(num_param) - f).is_nan());
  }

  res = f3 - f3;
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK(res.is_nan());
  } else {
    BOOST_CHECK_EQUAL(res(0,0), 0);
    BOOST_CHECK_EQUAL(res(0,1), 0);
    BOOST_CHECK_EQUAL(res(0,2), 0);
  }
  res = f3 - f4;
  BOOST_CHECK_EQUAL(res(0,0), f3(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f3(0,1));
  BOOST_CHECK_EQUAL(res(0,2), f3(0,2));

  res = f + f2;
  BOOST_CHECK_EQUAL(res(0,0), -5);
  BOOST_CHECK_EQUAL(res(0,1), 2);
  BOOST_CHECK_EQUAL(res(0,2), 0);

  res = f + f3;
  BOOST_CHECK_EQUAL(res(0,0), f3(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f3(0,1));
  BOOST_CHECK_EQUAL(res(0,2), f3(0,2));

  res = T(5) + f;
  BOOST_CHECK_EQUAL(res(0,0), -5);
  BOOST_CHECK_EQUAL(res(0,1), 5);
  BOOST_CHECK_EQUAL(res(0,2), 6);

  res = f + T(5);
  BOOST_CHECK_EQUAL(res(0,0), -5);
  BOOST_CHECK_EQUAL(res(0,1), 5);
  BOOST_CHECK_EQUAL(res(0,2), 6);

  BOOST_CHECK((f + F::inf(num_param)).is_plus_inf());
  BOOST_CHECK((F::inf(num_param) + f).is_plus_inf());
  BOOST_CHECK((f + F::minus_inf(num_param)).is_minus_inf());
  BOOST_CHECK((F::minus_inf(num_param) + f).is_minus_inf());
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK((f + F::nan(num_param)).is_nan());
    BOOST_CHECK((F::nan(num_param) + f).is_nan());
  }

  res = f3 + f4;
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK(res.is_nan());
  } else {
    BOOST_CHECK_EQUAL(res(0,0), 0);
    BOOST_CHECK_EQUAL(res(0,1), 0);
    BOOST_CHECK_EQUAL(res(0,2), 0);
  }
  res = f3 + f3;
  BOOST_CHECK_EQUAL(res(0,0), f3(0,0));
  BOOST_CHECK_EQUAL(res(0,1), f3(0,1));
  BOOST_CHECK_EQUAL(res(0,2), f3(0,2));

  res = f * f2;
  BOOST_CHECK_EQUAL(res(0,0), -50);
  BOOST_CHECK_EQUAL(res(0,1), 0);
  BOOST_CHECK_EQUAL(res(0,2), -1);

  res = f * f3;
  BOOST_CHECK_EQUAL(res(0,0), f4(0,0));
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }
  BOOST_CHECK_EQUAL(res(0,2), f3(0,2));

  res = T(5) * f;
  BOOST_CHECK_EQUAL(res(0,0), -50);
  BOOST_CHECK_EQUAL(res(0,1), 0);
  BOOST_CHECK_EQUAL(res(0,2), 5);

  res = f * T(5);
  BOOST_CHECK_EQUAL(res(0,0), -50);
  BOOST_CHECK_EQUAL(res(0,1),  0);
  BOOST_CHECK_EQUAL(res(0,2), 5);

  res = f * F::inf(num_param);
  BOOST_CHECK_EQUAL(res(0,0), -F::T_inf);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }
  BOOST_CHECK_EQUAL(res(0,2), F::T_inf);

  res = F::inf(num_param) * f;
  BOOST_CHECK_EQUAL(res(0,0), -F::T_inf);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }
  BOOST_CHECK_EQUAL(res(0,2), F::T_inf);

  res = f * F::minus_inf(num_param);
  BOOST_CHECK_EQUAL(res(0,0), F::T_inf);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }
  BOOST_CHECK_EQUAL(res(0,2), -F::T_inf);

  res = F::minus_inf(num_param) * f;
  BOOST_CHECK_EQUAL(res(0,0), F::T_inf);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }
  BOOST_CHECK_EQUAL(res(0,2), -F::T_inf);

  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    res = f * F::nan(num_param);
    BOOST_CHECK(res.is_nan());
    res = F::nan(num_param) * f;
    BOOST_CHECK(res.is_nan());
  }

  res = f3 * f3;
  BOOST_CHECK(res.is_plus_inf());
  res = f3 * f4;
  BOOST_CHECK(res.is_minus_inf());

  res = f / f2;
  BOOST_CHECK_EQUAL(res(0,0), -2);
  BOOST_CHECK_EQUAL(res(0,1), 0);
  BOOST_CHECK_EQUAL(res(0,2), -1);

  res = f / f3;
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 0);
  BOOST_CHECK_EQUAL(res(0,2), 0);

  res = f3 / f;
  BOOST_CHECK_EQUAL(res(0,0), f4(0,0));
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }
  BOOST_CHECK_EQUAL(res(0,2), f3(0,2));

  res = T(5) / f;
  BOOST_CHECK_EQUAL(res(0,0), static_cast<T>(-0.5));
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }
  BOOST_CHECK_EQUAL(res(0,2), 5);

  res = f / T(5);
  BOOST_CHECK_EQUAL(res(0,0), -2);
  BOOST_CHECK_EQUAL(res(0,1),  0);
  BOOST_CHECK_EQUAL(res(0,2), static_cast<T>(1) / static_cast<T>(5)); //to avoid precision error

  res = f / F::inf(num_param);
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 0);
  BOOST_CHECK_EQUAL(res(0,2), 0);
  res = F::inf(num_param) / f;
  BOOST_CHECK_EQUAL(res(0,0), -F::T_inf);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }
  BOOST_CHECK_EQUAL(res(0,2), F::T_inf);

  res = f / F::minus_inf(num_param);
  BOOST_CHECK_EQUAL(res(0,0), 0);
  BOOST_CHECK_EQUAL(res(0,1), 0);
  BOOST_CHECK_EQUAL(res(0,2), 0);
  res = F::minus_inf(num_param) / f;
  BOOST_CHECK_EQUAL(res(0,0), F::T_inf);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    BOOST_CHECK(std::isnan(res(0,1)));
  } else {
    BOOST_CHECK_EQUAL(res(0,1), 0);
  }
  BOOST_CHECK_EQUAL(res(0,2), -F::T_inf);

  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    res = f / F::nan(num_param);
    BOOST_CHECK(res.is_nan());
    res = F::nan(num_param) / f;
    BOOST_CHECK(res.is_nan());
  }

  res = f3 / f3;
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK(res.is_nan());
  } else {
    BOOST_CHECK_EQUAL(res(0,0), 0);
    BOOST_CHECK_EQUAL(res(0,1), 0);
    BOOST_CHECK_EQUAL(res(0,2), 0);
  }
  res = f3 / f4;
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK(res.is_nan());
  } else {
    BOOST_CHECK_EQUAL(res(0,0), 0);
    BOOST_CHECK_EQUAL(res(0,1), 0);
    BOOST_CHECK_EQUAL(res(0,2), 0);
  }
  res = f / F({0, 0, 0});
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    BOOST_CHECK(res.is_nan());
  } else {
    BOOST_CHECK_EQUAL(res(0,0), 0);
    BOOST_CHECK_EQUAL(res(0,1), 0);
    BOOST_CHECK_EQUAL(res(0,2), 0);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_operators, T, list_of_tested_variants)
{
  test_operators<Multi_parameter_filtration<T>, T>();
  test_operators<Multi_parameter_filtration<T, false, true>, T>();

  test_operators<Dynamic_multi_parameter_filtration<T>, T>();
  test_operators<Dynamic_multi_parameter_filtration<T, false, true>, T>();
}

template <class F, typename T>
void test_modifiers(){
  const int num_param = 3;
  std::vector<T> v;

  if constexpr (F::ensures_1_criticality()) {
    v = {0, 1, 2};
  } else {
    v = {-3, 1, 7, 0, 1, 2};
  }

  F f(v.begin(), v.end(), num_param);
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 1);
    BOOST_CHECK_EQUAL(f(0, 2), 2);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), -3);
    BOOST_CHECK_EQUAL(f(0, 1), 1);
    BOOST_CHECK_EQUAL(f(0, 2), 7);
    BOOST_CHECK_EQUAL(f(1, 0), 0);
    BOOST_CHECK_EQUAL(f(1, 1), 1);
    BOOST_CHECK_EQUAL(f(1, 2), 2);
  }

  f.push_to_least_common_upper_bound({-1, 5, 6});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 5);
    BOOST_CHECK_EQUAL(f(0, 2), 6);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), -1);
    BOOST_CHECK_EQUAL(f(0, 1), 5);
    BOOST_CHECK_EQUAL(f(0, 2), 7);
    BOOST_CHECK_EQUAL(f(1, 0), 0);
    BOOST_CHECK_EQUAL(f(1, 1), 5);
    BOOST_CHECK_EQUAL(f(1, 2), 6);
  }

  f.push_to_least_common_upper_bound({-1, -5, -6});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 5);
    BOOST_CHECK_EQUAL(f(0, 2), 6);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), -1);
    BOOST_CHECK_EQUAL(f(0, 1), 5);
    BOOST_CHECK_EQUAL(f(0, 2), 7);
    BOOST_CHECK_EQUAL(f(1, 0), 0);
    BOOST_CHECK_EQUAL(f(1, 1), 5);
    BOOST_CHECK_EQUAL(f(1, 2), 6);
  }

  f.push_to_least_common_upper_bound(F::minus_inf(num_param));
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f(0, 0), 0);
    BOOST_CHECK_EQUAL(f(0, 1), 5);
    BOOST_CHECK_EQUAL(f(0, 2), 6);
  } else {
    BOOST_CHECK_EQUAL(f(0, 0), -1);
    BOOST_CHECK_EQUAL(f(0, 1), 5);
    BOOST_CHECK_EQUAL(f(0, 2), 7);
    BOOST_CHECK_EQUAL(f(1, 0), 0);
    BOOST_CHECK_EQUAL(f(1, 1), 5);
    BOOST_CHECK_EQUAL(f(1, 2), 6);
  }

  f.push_to_least_common_upper_bound(F::inf(num_param));
  BOOST_CHECK(f.is_plus_inf());

  if constexpr (std::numeric_limits<F>::has_quiet_NaN){
    f.push_to_least_common_upper_bound(F::nan(num_param));
    BOOST_CHECK(f.is_plus_inf());
  }

  f.pull_to_greatest_common_lower_bound({-1, 5, 6});
  BOOST_CHECK_EQUAL(f(0, 0), -1);
  BOOST_CHECK_EQUAL(f(0, 1), 5);
  BOOST_CHECK_EQUAL(f(0, 2), 6);

  f.pull_to_greatest_common_lower_bound({1, 8, 9});
  BOOST_CHECK_EQUAL(f(0, 0), -1);
  BOOST_CHECK_EQUAL(f(0, 1), 5);
  BOOST_CHECK_EQUAL(f(0, 2), 6);

  f.pull_to_greatest_common_lower_bound(F::inf(num_param));
  BOOST_CHECK_EQUAL(f(0, 0), -1);
  BOOST_CHECK_EQUAL(f(0, 1), 5);
  BOOST_CHECK_EQUAL(f(0, 2), 6);

  f.pull_to_greatest_common_lower_bound(F::minus_inf(num_param));
  BOOST_CHECK(f.is_minus_inf());

  if constexpr (std::numeric_limits<F>::has_quiet_NaN){
    f.pull_to_greatest_common_lower_bound(F::nan(num_param));
    BOOST_CHECK(f.is_minus_inf());
  }

  std::vector<std::vector<int> > grid = {{0, 2, 4, 8}, {0, 3, 6, 9}, {0, 4, 8, 16}};

  f.push_to_least_common_upper_bound({1, 7, 5});
  f.project_onto_grid(grid, true);
  BOOST_CHECK_EQUAL(f(0, 0), 1);
  BOOST_CHECK_EQUAL(f(0, 1), 3);
  BOOST_CHECK_EQUAL(f(0, 2), 2);

  f.push_to_least_common_upper_bound({1, 7, 5});
  f.project_onto_grid(grid, false);
  BOOST_CHECK_EQUAL(f(0, 0), 2);
  BOOST_CHECK_EQUAL(f(0, 1), 9);
  BOOST_CHECK_EQUAL(f(0, 2), 8);

  if constexpr (!F::ensures_1_criticality()) {
    f.set_num_generators(3);
    BOOST_CHECK_EQUAL(f.num_parameters(), 3);
    BOOST_CHECK_EQUAL(f.num_generators(), 3);
    BOOST_CHECK_EQUAL(f.num_entries(), 9);
    BOOST_CHECK_EQUAL(f(0, 0), 2);
    BOOST_CHECK_EQUAL(f(0, 1), 9);
    BOOST_CHECK_EQUAL(f(0, 2), 8);
    BOOST_CHECK_EQUAL(f(1, 0), -F::T_inf);
    BOOST_CHECK_EQUAL(f(1, 1), -F::T_inf);
    BOOST_CHECK_EQUAL(f(1, 2), -F::T_inf);
    BOOST_CHECK_EQUAL(f(2, 0), -F::T_inf);
    BOOST_CHECK_EQUAL(f(2, 1), -F::T_inf);
    BOOST_CHECK_EQUAL(f(2, 2), -F::T_inf);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_modifiers, T, list_of_tested_variants)
{
  test_modifiers<Multi_parameter_filtration<T>, T>();
  test_modifiers<Multi_parameter_filtration<T, false, true>, T>();

  test_modifiers<Dynamic_multi_parameter_filtration<T>, T>();
  test_modifiers<Dynamic_multi_parameter_filtration<T, false, true>, T>();
}

template <class F, typename T>
void test_add_generators(){
  const int num_param = 3;

  F f({0, 1, 2});
  BOOST_CHECK_EQUAL(f.num_generators(), 1);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), 0);
  BOOST_CHECK_EQUAL(f(0, 1), 1);
  BOOST_CHECK_EQUAL(f(0, 2), 2);

  bool res = f.add_generator({-3, 1, 7});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), -3);
  BOOST_CHECK_EQUAL(f(0, 1), 1);
  BOOST_CHECK_EQUAL(f(0, 2), 7);
  BOOST_CHECK_EQUAL(f(1, 0), 0);
  BOOST_CHECK_EQUAL(f(1, 1), 1);
  BOOST_CHECK_EQUAL(f(1, 2), 2);

  res = f.add_generator({-1, -2, -3});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), -3);
  BOOST_CHECK_EQUAL(f(0, 1), 1);
  BOOST_CHECK_EQUAL(f(0, 2), 7);
  BOOST_CHECK_EQUAL(f(1, 0), -1);
  BOOST_CHECK_EQUAL(f(1, 1), -2);
  BOOST_CHECK_EQUAL(f(1, 2), -3);

  res = f.add_generator({8, 9, 10});
  BOOST_CHECK(!res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), -3);
  BOOST_CHECK_EQUAL(f(0, 1), 1);
  BOOST_CHECK_EQUAL(f(0, 2), 7);
  BOOST_CHECK_EQUAL(f(1, 0), -1);
  BOOST_CHECK_EQUAL(f(1, 1), -2);
  BOOST_CHECK_EQUAL(f(1, 2), -3);

  res = f.add_generator(std::vector<T>(num_param, F::T_inf));
  BOOST_CHECK(!res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), -3);
  BOOST_CHECK_EQUAL(f(0, 1), 1);
  BOOST_CHECK_EQUAL(f(0, 2), 7);
  BOOST_CHECK_EQUAL(f(1, 0), -1);
  BOOST_CHECK_EQUAL(f(1, 1), -2);
  BOOST_CHECK_EQUAL(f(1, 2), -3);

  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    res = f.add_generator(std::vector<T>(num_param, std::numeric_limits<T>::quiet_NaN()));
    BOOST_CHECK(!res);
    BOOST_CHECK_EQUAL(f.num_generators(), 2);
    BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
    BOOST_CHECK_EQUAL(f(0, 0), -3);
    BOOST_CHECK_EQUAL(f(0, 1), 1);
    BOOST_CHECK_EQUAL(f(0, 2), 7);
    BOOST_CHECK_EQUAL(f(1, 0), -1);
    BOOST_CHECK_EQUAL(f(1, 1), -2);
    BOOST_CHECK_EQUAL(f(1, 2), -3);
  }

  res = f.add_generator(std::vector<T>(num_param, -F::T_inf));
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 1);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), -F::T_inf);

  std::vector<T> v{
      0, 1, 2,
      F::T_inf, F::T_inf, F::T_inf,
      0, 1, 2,
      std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN(),
      -F::T_inf, -F::T_inf, -F::T_inf};

  F f2(v.begin(), v.end(), num_param);
  f2.remove_empty_generators(false);
  BOOST_CHECK_EQUAL(f2.num_generators(), 5);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(f2(0, 0), -F::T_inf);
    BOOST_CHECK_EQUAL(f2(0, 1), -F::T_inf);
    BOOST_CHECK_EQUAL(f2(0, 2), -F::T_inf);
    BOOST_CHECK_EQUAL(f2(1, 0), 0);
    BOOST_CHECK_EQUAL(f2(1, 1), 1);
    BOOST_CHECK_EQUAL(f2(1, 2), 2);
    BOOST_CHECK_EQUAL(f2(2, 0), 0);
    BOOST_CHECK_EQUAL(f2(2, 1), 1);
    BOOST_CHECK_EQUAL(f2(2, 2), 2);
    BOOST_CHECK_EQUAL(f2(3, 0), F::T_inf);
    BOOST_CHECK_EQUAL(f2(3, 1), F::T_inf);
    BOOST_CHECK_EQUAL(f2(3, 2), F::T_inf);
    BOOST_CHECK(std::isnan(f2(4, 0)));
    BOOST_CHECK(std::isnan(f2(4, 1)));
    BOOST_CHECK(std::isnan(f2(4, 2)));
  } else {
    BOOST_CHECK_EQUAL(f2(0, 0), -F::T_inf);
    BOOST_CHECK_EQUAL(f2(0, 1), -F::T_inf);
    BOOST_CHECK_EQUAL(f2(0, 2), -F::T_inf);
    BOOST_CHECK_EQUAL(f2(1, 0), 0);
    BOOST_CHECK_EQUAL(f2(1, 1), 0);
    BOOST_CHECK_EQUAL(f2(1, 2), 0);
    BOOST_CHECK_EQUAL(f2(2, 0), 0);
    BOOST_CHECK_EQUAL(f2(2, 1), 1);
    BOOST_CHECK_EQUAL(f2(2, 2), 2);
    BOOST_CHECK_EQUAL(f2(3, 0), 0);
    BOOST_CHECK_EQUAL(f2(3, 1), 1);
    BOOST_CHECK_EQUAL(f2(3, 2), 2);
    BOOST_CHECK_EQUAL(f2(4, 0), F::T_inf);
    BOOST_CHECK_EQUAL(f2(4, 1), F::T_inf);
    BOOST_CHECK_EQUAL(f2(4, 2), F::T_inf);
  }

  F f3(v.begin(), v.end(), num_param);
  f3.remove_empty_generators(true);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(f3.num_generators(), 2);
    BOOST_CHECK_EQUAL(f3(0, 0), 0);
    BOOST_CHECK_EQUAL(f3(0, 1), 1);
    BOOST_CHECK_EQUAL(f3(0, 2), 2);
    BOOST_CHECK_EQUAL(f3(1, 0), 0);
    BOOST_CHECK_EQUAL(f3(1, 1), 1);
    BOOST_CHECK_EQUAL(f3(1, 2), 2);
  } else {
    BOOST_CHECK_EQUAL(f3.num_generators(), 3);
    BOOST_CHECK_EQUAL(f3(0, 0), 0);
    BOOST_CHECK_EQUAL(f3(0, 1), 0);
    BOOST_CHECK_EQUAL(f3(0, 2), 0);
    BOOST_CHECK_EQUAL(f3(1, 0), 0);
    BOOST_CHECK_EQUAL(f3(1, 1), 1);
    BOOST_CHECK_EQUAL(f3(1, 2), 2);
    BOOST_CHECK_EQUAL(f3(2, 0), 0);
    BOOST_CHECK_EQUAL(f3(2, 1), 1);
    BOOST_CHECK_EQUAL(f3(2, 2), 2);
  }

  F f4(v.begin(), v.end(), num_param);
  f4.simplify();
  BOOST_CHECK_EQUAL(f4.num_generators(), 1);
  BOOST_CHECK(f4.is_minus_inf());

  v.pop_back();
  F f5(v.begin(), v.end(), num_param);
  f5.simplify();
  BOOST_CHECK_EQUAL(f5.num_generators(), 1);
  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK_EQUAL(f5(0, 0), 0);
    BOOST_CHECK_EQUAL(f5(0, 1), 1);
    BOOST_CHECK_EQUAL(f5(0, 2), 2);
  } else {
    BOOST_CHECK_EQUAL(f5(0, 0), 0);
    BOOST_CHECK_EQUAL(f5(0, 1), 0);
    BOOST_CHECK_EQUAL(f5(0, 2), 0);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_add_generators, T, list_of_tested_variants)
{
  test_add_generators<Multi_parameter_filtration<T>, T>();

  Multi_parameter_filtration<T, false, true> f({0, 1, 2});
  BOOST_CHECK_THROW(f.add_generator({-3, 1, 7}), std::logic_error);

  test_add_generators<Dynamic_multi_parameter_filtration<T>, T>();

  Dynamic_multi_parameter_filtration<T, false, true> f2({0, 1, 2});
  BOOST_CHECK_THROW(f2.add_generator({-3, 1, 7}), std::logic_error);
}

template <class F, typename T>
void test_friends(){
  F f({0, 1, 2});

  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(5))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, std::initializer_list<T>{2, 3, 5}),
                    static_cast<T>(std::sqrt(T(17))));
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {2, 3, 5, 9}), 13);
  BOOST_CHECK(factorize_below(f) == f);
  BOOST_CHECK(factorize_above(f) == f);

  if constexpr (!F::ensures_1_criticality()) {
    f.add_guaranteed_generator({2, 0, 4});
    BOOST_CHECK_EQUAL(f.num_generators(), 2);
    BOOST_CHECK_EQUAL(f.num_parameters(), 3);

    BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(25))));
    BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, std::initializer_list<T>{2, 3, 5}),
                      static_cast<T>(std::sqrt(T(10))));
    BOOST_CHECK_EQUAL(compute_linear_projection(f, {2, 3, 5, 9}), 13);
    BOOST_CHECK(factorize_below(f) == F({0, 0, 2}));
    BOOST_CHECK(factorize_above(f) == F({2, 1, 4}));

    if constexpr (std::numeric_limits<T>::has_quiet_NaN){
      T nan = std::numeric_limits<T>::quiet_NaN();
      std::vector<T> v = {0, nan, 2, 2, nan, 4};
      F f2(v.begin(), v.end(), 3);

      BOOST_CHECK(std::isnan(compute_norm(f2)));
      BOOST_CHECK(std::isnan(compute_euclidean_distance_to(f2, std::initializer_list<T>{2, 3, 5})));
      BOOST_CHECK(std::isnan(compute_linear_projection(f2, {2, 3, 5, 9})));
      auto bf2 = factorize_below(f2);
      BOOST_CHECK_EQUAL(bf2(0, 0), 0);
      BOOST_CHECK(std::isnan(bf2(0, 1)));
      BOOST_CHECK_EQUAL(bf2(0, 2), 2);
      auto af2 = factorize_above(f2);
      BOOST_CHECK_EQUAL(af2(0, 0), 2);
      BOOST_CHECK(std::isnan(af2(0, 1)));
      BOOST_CHECK_EQUAL(af2(0, 2), 4);
    }
  }

  f(0,0) = 1;
  f(0,1) = 7;
  f(0,2) = 5;

  std::vector<std::vector<int> > grid = {{0, 2, 4, 8}, {0, 3, 6, 9}, {0, 4, 8, 16}};
  auto res = compute_coordinates_in_grid(f, grid);
  BOOST_CHECK_EQUAL(res.num_parameters(), 3);
  BOOST_CHECK_EQUAL(f.num_parameters(), 3);
  BOOST_CHECK_EQUAL(res(0, 0), 1);
  BOOST_CHECK_EQUAL(res(0, 1), 3);
  BOOST_CHECK_EQUAL(res(0, 2), 2);

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(res.num_generators(), 1);
    BOOST_CHECK_EQUAL(f.num_generators(), 1);
  } else {
    BOOST_CHECK_EQUAL(res.num_generators(), 2);
    BOOST_CHECK_EQUAL(f.num_generators(), 2);
    BOOST_CHECK_EQUAL(res(1, 0), 1);
    BOOST_CHECK_EQUAL(res(1, 1), 0);
    BOOST_CHECK_EQUAL(res(1, 2), 1);
  }

  res = evaluate_coordinates_in_grid(res, grid);
  BOOST_CHECK_EQUAL(res.num_generators(), 1);
  BOOST_CHECK_EQUAL(res.num_parameters(), 3);

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(res(0, 0), 2);
    BOOST_CHECK_EQUAL(res(0, 1), 9);
    BOOST_CHECK_EQUAL(res(0, 2), 8);
  } else {
    BOOST_CHECK_EQUAL(res(0, 0), 2);
    BOOST_CHECK_EQUAL(res(0, 1), 0);
    BOOST_CHECK_EQUAL(res(0, 2), 4);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_friends, T, list_of_tested_variants)
{
  test_friends<Multi_parameter_filtration<T>, T>();
  test_friends<Multi_parameter_filtration<T, false, true>, T>();

  test_friends<Dynamic_multi_parameter_filtration<T>, T>();
  test_friends<Dynamic_multi_parameter_filtration<T, false, true>, T>();
}

template <class F, typename T>
void test_unify_intersect(){
  const int num_param = 2;

  std::vector<T> v1 = {0,5,2,3,5,2};
  F f1(v1.begin(), v1.end(), num_param);

  std::vector<T> v2 = {1,4,4,1};
  F f2(v2.begin(), v2.end(), num_param);

  bool modified = unify_lifetimes(f1, f2);
  BOOST_CHECK(modified);
  BOOST_CHECK(f1.num_parameters() == num_param);
  BOOST_CHECK(f1.num_generators() == 4);
  BOOST_CHECK_EQUAL(f1(0,0), 0);
  BOOST_CHECK_EQUAL(f1(0,1), 5);
  BOOST_CHECK_EQUAL(f1(1,0), 1);
  BOOST_CHECK_EQUAL(f1(1,1), 4);
  BOOST_CHECK_EQUAL(f1(2,0), 2);
  BOOST_CHECK_EQUAL(f1(2,1), 3);
  BOOST_CHECK_EQUAL(f1(3,0), 4);
  BOOST_CHECK_EQUAL(f1(3,1), 1);

  std::vector<T> v3 = {0,5,2,3,5,2};
  F f3(v3.begin(), v3.end(), num_param);

  modified = intersect_lifetimes(f3, f2);
  BOOST_CHECK(modified);
  BOOST_CHECK(f3.num_parameters() == num_param);
  BOOST_CHECK(f3.num_generators() == 4);
  BOOST_CHECK_EQUAL(f3(0,0), 1);
  BOOST_CHECK_EQUAL(f3(0,1), 5);
  BOOST_CHECK_EQUAL(f3(1,0), 2);
  BOOST_CHECK_EQUAL(f3(1,1), 4);
  BOOST_CHECK_EQUAL(f3(2,0), 4);
  BOOST_CHECK_EQUAL(f3(2,1), 3);
  BOOST_CHECK_EQUAL(f3(3,0), 5);
  BOOST_CHECK_EQUAL(f3(3,1), 2);

  modified = unify_lifetimes(f1, f3);
  BOOST_CHECK(!modified);
}

template <class F, typename T>
void test_unify_intersect_1_critical(){
  const int num_param = 3;

  std::vector<T> v1 = {0,5,2};
  F f1(v1.begin(), v1.end(), num_param);

  std::vector<T> v2 = {1,8,4};
  F f2(v2.begin(), v2.end(), num_param);

  bool modified = unify_lifetimes(f1, f2);
  BOOST_CHECK(!modified);
  BOOST_CHECK(f1.num_parameters() == num_param);
  BOOST_CHECK(f1.num_generators() == 1);
  BOOST_CHECK_EQUAL(f1(0,0), 0);
  BOOST_CHECK_EQUAL(f1(0,1), 5);
  BOOST_CHECK_EQUAL(f1(0,2), 2);

  modified = intersect_lifetimes(f1, f2);
  BOOST_CHECK(modified);
  BOOST_CHECK(f1.num_parameters() == num_param);
  BOOST_CHECK(f1.num_generators() == 1);
  BOOST_CHECK_EQUAL(f1(0,0), 1);
  BOOST_CHECK_EQUAL(f1(0,1), 8);
  BOOST_CHECK_EQUAL(f1(0,2), 4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_unify_intersect, T, list_of_tested_variants)
{
  test_unify_intersect<Multi_parameter_filtration<T>, T>();
  test_unify_intersect_1_critical<Multi_parameter_filtration<T, false, true>, T>();

  test_unify_intersect<Dynamic_multi_parameter_filtration<T>, T>();
  test_unify_intersect_1_critical<Dynamic_multi_parameter_filtration<T, false, true>, T>();
}

template <class F, typename T>
void test_serialize(){
  std::vector<T> v = {0,5,2,3,5,2};
  F f1(v.begin(), v.end(), 2);
  BOOST_CHECK(f1.num_parameters() == 2);
  BOOST_CHECK(f1.num_generators() == 3);
  BOOST_CHECK_EQUAL(f1(0,0), 0);
  BOOST_CHECK_EQUAL(f1(0,1), 5);
  BOOST_CHECK_EQUAL(f1(1,0), 2);
  BOOST_CHECK_EQUAL(f1(1,1), 3);
  BOOST_CHECK_EQUAL(f1(2,0), 5);
  BOOST_CHECK_EQUAL(f1(2,1), 2);
  F f2(v.begin(), v.end(), 3);
  BOOST_CHECK(f2.num_parameters() == 3);
  BOOST_CHECK(f2.num_generators() == 2);
  BOOST_CHECK_EQUAL(f2(0,0), 0);
  BOOST_CHECK_EQUAL(f2(0,1), 5);
  BOOST_CHECK_EQUAL(f2(0,2), 2);
  BOOST_CHECK_EQUAL(f2(1,0), 3);
  BOOST_CHECK_EQUAL(f2(1,1), 5);
  BOOST_CHECK_EQUAL(f2(1,2), 2);

  char* buffer = new char[256];
  std::size_t serializationSize = get_serialization_size_of(f1);
  
  char* ptr = buffer;
  ptr = serialize_value_to_char_buffer(f1, ptr);
  BOOST_CHECK_EQUAL((void*)ptr, (void*)(buffer + serializationSize));

  const char* c_ptr = buffer;
  F f3(0);
  c_ptr = deserialize_value_from_char_buffer(f3, c_ptr);
  BOOST_CHECK_EQUAL((void*)c_ptr, (void*)(buffer + serializationSize));
  BOOST_CHECK(f3.num_parameters() == 2);
  BOOST_CHECK(f3.num_generators() == 3);
  BOOST_CHECK_EQUAL(f3(0,0), 0);
  BOOST_CHECK_EQUAL(f3(0,1), 5);
  BOOST_CHECK_EQUAL(f3(1,0), 2);
  BOOST_CHECK_EQUAL(f3(1,1), 3);
  BOOST_CHECK_EQUAL(f3(2,0), 5);
  BOOST_CHECK_EQUAL(f3(2,1), 2);

  serializationSize = get_serialization_size_of(f2);
  
  ptr = buffer;
  ptr = serialize_value_to_char_buffer(f2, ptr);
  BOOST_CHECK_EQUAL((void*)ptr, (void*)(buffer + serializationSize));

  c_ptr = buffer;
  F f4(0);
  c_ptr = deserialize_value_from_char_buffer(f4, c_ptr);
  BOOST_CHECK_EQUAL((void*)c_ptr, (void*)(buffer + serializationSize));
  BOOST_CHECK(f4.num_parameters() == 3);
  BOOST_CHECK(f4.num_generators() == 2);
  BOOST_CHECK_EQUAL(f4(0,0), 0);
  BOOST_CHECK_EQUAL(f4(0,1), 5);
  BOOST_CHECK_EQUAL(f4(0,2), 2);
  BOOST_CHECK_EQUAL(f4(1,0), 3);
  BOOST_CHECK_EQUAL(f4(1,1), 5);
  BOOST_CHECK_EQUAL(f4(1,2), 2);

  delete[] buffer;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_serialize, T, list_of_tested_variants)
{
  test_serialize<Multi_parameter_filtration<T>, T>();

  test_serialize<Dynamic_multi_parameter_filtration<T>, T>();
}

template <class F, typename T>
void test_co(){
  F f;
  BOOST_CHECK(f.num_parameters() == 2);
  BOOST_CHECK(f.num_generators() == 1);
  BOOST_CHECK_EQUAL(f(0,0), F::T_inf);
  BOOST_CHECK_EQUAL(f(0,1), F::T_inf);

  BOOST_CHECK(f.is_plus_inf());
  BOOST_CHECK(!f.is_minus_inf());
  BOOST_CHECK(!f.is_nan());
  BOOST_CHECK(!f.is_finite());

  F f6 = F::minus_inf(3);
  bool change = f6.add_generator(std::vector<T>(3, F::T_inf));
  BOOST_CHECK(change);
  BOOST_CHECK(f6.is_plus_inf());

  if constexpr (F::ensures_1_criticality()) {
    std::vector<T> v = {0, 1, 2};
    F f2(v.begin(), v.end(), 3);
    BOOST_CHECK_EQUAL(compute_linear_projection(f2, {2,3,5,9}), 13);
  } else {
    std::vector<T> v = {0, 1, 2, 2, 0, 4};
    F f2(v.begin(), v.end(), 3);
    BOOST_CHECK_EQUAL(compute_linear_projection(f2, {2,3,5,9}), 24);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_co, T, list_of_tested_variants)
{
  test_co<Multi_parameter_filtration<T, true>, T>();
  test_co<Multi_parameter_filtration<T, true, true>, T>();

  test_co<Dynamic_multi_parameter_filtration<T, true>, T>();
  test_co<Dynamic_multi_parameter_filtration<T, true, true>, T>();
}

template <class F, typename T>
void test_numerical_limits(){
  const int num_param = 3;

  BOOST_CHECK(std::numeric_limits<F>::has_infinity);

  BOOST_CHECK_THROW(std::numeric_limits<F>::infinity(), std::logic_error);
  BOOST_CHECK_THROW(std::numeric_limits<F>::minus_infinity(), std::logic_error);
  BOOST_CHECK_THROW(std::numeric_limits<F>::quiet_NaN(), std::logic_error);
  BOOST_CHECK_THROW(std::numeric_limits<F>::max(), std::logic_error);

  BOOST_CHECK(std::numeric_limits<F>::infinity(num_param).is_plus_inf());
  BOOST_CHECK(std::numeric_limits<F>::minus_infinity(num_param).is_minus_inf());
  
  auto max = std::numeric_limits<F>::max(num_param);
  BOOST_CHECK_EQUAL(max(0,0), std::numeric_limits<T>::max());
  BOOST_CHECK_EQUAL(max(0,1), std::numeric_limits<T>::max());
  BOOST_CHECK_EQUAL(max(0,2), std::numeric_limits<T>::max());

  if constexpr (std::numeric_limits<T>::has_quiet_NaN){
    BOOST_CHECK(std::numeric_limits<F>::has_quiet_NaN);
    BOOST_CHECK(std::numeric_limits<F>::quiet_NaN(num_param).is_nan());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_numerical_limits, T, list_of_tested_variants)
{
  test_numerical_limits<Multi_parameter_filtration<T>, T>();
  test_numerical_limits<Multi_parameter_filtration<T, false, true>, T>();

  test_numerical_limits<Dynamic_multi_parameter_filtration<T>, T>();
  test_numerical_limits<Dynamic_multi_parameter_filtration<T, false, true>, T>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_converters, T, list_of_tested_variants)
{
  std::vector<T> v = {0, 5, 1, 4, 2, 3, 3, 2};
  Dynamic_multi_parameter_filtration<T> f0(v.begin(), v.end(), 2);
  BOOST_CHECK(f0.num_parameters() == 2);
  BOOST_CHECK(f0.num_generators() == 4);
  BOOST_CHECK_EQUAL(f0(0,0), 0);
  BOOST_CHECK_EQUAL(f0(0,1), 5);
  BOOST_CHECK_EQUAL(f0(1,0), 1);
  BOOST_CHECK_EQUAL(f0(1,1), 4);
  BOOST_CHECK_EQUAL(f0(2,0), 2);
  BOOST_CHECK_EQUAL(f0(2,1), 3);
  BOOST_CHECK_EQUAL(f0(3,0), 3);
  BOOST_CHECK_EQUAL(f0(3,1), 2);

  Multi_parameter_filtration<T> f1 = f0.convert_to_multi_parameter_filtration();
  BOOST_CHECK(f1.num_parameters() == 2);
  BOOST_CHECK(f1.num_generators() == 4);
  BOOST_CHECK_EQUAL(f1(0,0), 0);
  BOOST_CHECK_EQUAL(f1(0,1), 5);
  BOOST_CHECK_EQUAL(f1(1,0), 1);
  BOOST_CHECK_EQUAL(f1(1,1), 4);
  BOOST_CHECK_EQUAL(f1(2,0), 2);
  BOOST_CHECK_EQUAL(f1(2,1), 3);
  BOOST_CHECK_EQUAL(f1(3,0), 3);
  BOOST_CHECK_EQUAL(f1(3,1), 2);
}

