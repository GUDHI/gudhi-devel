/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstddef>    //std::size_t
#include <limits>     //std::numerical_limits
#include <stdexcept>  //std::logic_error, std::out_of_range
#include <utility>    //std::swap, std::move
#include <vector>
#include <initializer_list>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_filtration"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_filtration/multi_filtration_utils.h>  // _is_nan with Windows fix
#include <gudhi/Multi_filtration/Degree_bifiltration.h>
#include <gudhi/Multi_filtration/Flat_array_filtration.h>
#include <gudhi/Multi_filtration/Nested_array_filtration.h>
#include <gudhi/Multi_filtration/multi_filtration_conversions.h>
#include <gudhi/Multi_filtration/multi_filtration_products.h>
#include <gudhi/Multi_parameter_filtration_value.h>

using namespace Gudhi::multi_filtration;

// declaration needed pre C++20
template <typename U, class MultiFiltrationValue, class CoefficientRange>
U compute_linear_projection();
template <typename U, class MultiFiltrationValue>
U compute_euclidean_distance_to();
template <typename U, class MultiFiltrationValue>
U compute_norm();

using list_of_tested_variants = boost::mpl::list<double, float, int>;

template <class F, typename T>
void test_constructors() {
  const int numParam = 2;

  F f0;
  BOOST_CHECK_EQUAL(f0.num_entries(), 2);
  BOOST_CHECK_EQUAL(f0.num_generators(), 1);
  BOOST_CHECK_EQUAL(f0.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f0(0, 1), 0);
  BOOST_CHECK_EQUAL(f0(0, 0), F::T_m_inf);

  F f1(numParam);
  BOOST_CHECK_EQUAL(f1.num_entries(), 2);
  BOOST_CHECK_EQUAL(f1.num_generators(), 1);
  BOOST_CHECK_EQUAL(f1.num_parameters(), 2);
  BOOST_CHECK_EQUAL((f1[{0, 1}]), 0);
  BOOST_CHECK_EQUAL((f1[{0, 0}]), F::T_m_inf);

  F f2(numParam, 0);
  BOOST_CHECK_EQUAL(f2.num_entries(), 2);
  BOOST_CHECK_EQUAL(f2.num_generators(), 1);
  BOOST_CHECK_EQUAL(f2.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f2(0, 1), 0);
  BOOST_CHECK_EQUAL(f2(0, 0), 0);

  F f3({1, 0, 2});
  BOOST_CHECK_EQUAL(f3.num_entries(), 2);
  BOOST_CHECK_EQUAL(f3.num_generators(), 1);
  BOOST_CHECK_EQUAL(f3.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f3(0, 1), 0);
  BOOST_CHECK_EQUAL(f3(0, 0), 1);

  std::vector<T> v{1, 0, 3, 1, 5, 2};
  F f4(v.begin(), v.end());
  BOOST_CHECK_EQUAL(f4.num_entries(), 2);
  BOOST_CHECK_EQUAL(f4.num_generators(), 1);
  BOOST_CHECK_EQUAL(f4.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f4(0, 1), 0);
  BOOST_CHECK_EQUAL(f4(0, 0), 1);

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_THROW(F f9(v.begin(), v.end(), numParam), std::logic_error);
    if constexpr (std::is_same_v<std::vector<T>, typename F::Underlying_container>) {
      BOOST_CHECK_THROW(F f5(v, 3), std::logic_error);
      BOOST_CHECK_THROW(F f6(std::move(v), 3), std::logic_error);
    }
  } else {
    F f9(v.begin(), v.end(), numParam);
    BOOST_CHECK_EQUAL(f9.num_entries(), 6);
    BOOST_CHECK_EQUAL(f9.num_generators(), 3);
    BOOST_CHECK_EQUAL(f9.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f9(0, 0), 1);
    BOOST_CHECK_EQUAL(f9(0, 1), 0);
    BOOST_CHECK_EQUAL(f9(1, 0), 3);
    BOOST_CHECK_EQUAL(f9(1, 1), 1);
    BOOST_CHECK_EQUAL(f9(2, 0), 5);
    BOOST_CHECK_EQUAL(f9(2, 1), 2);

    if constexpr (std::is_same_v<std::vector<T>, typename F::Underlying_container>) {
      F f5(v, numParam);
      BOOST_CHECK_EQUAL(f5.num_entries(), 12);
      BOOST_CHECK_EQUAL(f5.num_generators(), 6);
      BOOST_CHECK_EQUAL(f5.num_parameters(), 2);
      BOOST_CHECK_EQUAL((f5[{0, 0}]), 1);
      BOOST_CHECK_EQUAL((f5[{0, 1}]), 0);
      BOOST_CHECK_EQUAL((f5[{1, 0}]), 0);
      BOOST_CHECK_EQUAL((f5[{1, 1}]), 1);
      BOOST_CHECK_EQUAL((f5[{2, 0}]), 3);
      BOOST_CHECK_EQUAL((f5[{2, 1}]), 2);
      BOOST_CHECK_EQUAL((f5[{3, 0}]), 1);
      BOOST_CHECK_EQUAL((f5[{3, 1}]), 3);
      BOOST_CHECK_EQUAL((f5[{4, 0}]), 5);
      BOOST_CHECK_EQUAL((f5[{4, 1}]), 4);
      BOOST_CHECK_EQUAL((f5[{5, 0}]), 2);
      BOOST_CHECK_EQUAL((f5[{5, 1}]), 5);

      F f6(std::move(v), numParam);
      BOOST_CHECK(v.empty());
      BOOST_CHECK_EQUAL(f6.num_entries(), 12);
      BOOST_CHECK_EQUAL(f6.num_generators(), 6);
      BOOST_CHECK_EQUAL(f6.num_parameters(), 2);
      BOOST_CHECK_EQUAL(f6(0, 0), 1);
      BOOST_CHECK_EQUAL(f6(0, 1), 0);
      BOOST_CHECK_EQUAL(f6(1, 0), 0);
      BOOST_CHECK_EQUAL(f6(1, 1), 1);
      BOOST_CHECK_EQUAL(f6(2, 0), 3);
      BOOST_CHECK_EQUAL(f6(2, 1), 2);
      BOOST_CHECK_EQUAL(f6(3, 0), 1);
      BOOST_CHECK_EQUAL(f6(3, 1), 3);
      BOOST_CHECK_EQUAL(f6(4, 0), 5);
      BOOST_CHECK_EQUAL(f6(4, 1), 4);
      BOOST_CHECK_EQUAL(f6(5, 0), 2);
      BOOST_CHECK_EQUAL(f6(5, 1), 5);
    }

    F f7(f9);
    BOOST_CHECK_EQUAL(f7.num_entries(), 6);
    BOOST_CHECK_EQUAL(f7.num_generators(), 3);
    BOOST_CHECK_EQUAL(f7.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f7(0, 0), 1);
    BOOST_CHECK_EQUAL(f7(0, 1), 0);
    BOOST_CHECK_EQUAL(f7(1, 0), 3);
    BOOST_CHECK_EQUAL(f7(1, 1), 1);
    BOOST_CHECK_EQUAL(f7(2, 0), 5);
    BOOST_CHECK_EQUAL(f7(2, 1), 2);

    F f8(std::move(f9));
    BOOST_CHECK_EQUAL(f8.num_entries(), 6);
    BOOST_CHECK_EQUAL(f8.num_generators(), 3);
    BOOST_CHECK_EQUAL(f8.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f8(0, 0), 1);
    BOOST_CHECK_EQUAL(f8(0, 1), 0);
    BOOST_CHECK_EQUAL(f8(1, 0), 3);
    BOOST_CHECK_EQUAL(f8(1, 1), 1);
    BOOST_CHECK_EQUAL(f8(2, 0), 5);
    BOOST_CHECK_EQUAL(f8(2, 1), 2);

    swap(f0, f8);
    BOOST_CHECK_EQUAL(f8.num_entries(), 2);
    BOOST_CHECK_EQUAL(f8.num_generators(), 1);
    BOOST_CHECK_EQUAL(f8.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f8(0, 1), 0);
    BOOST_CHECK_EQUAL(f8(0, 0), F::T_m_inf);
    BOOST_CHECK_EQUAL(f0.num_entries(), 6);
    BOOST_CHECK_EQUAL(f0.num_generators(), 3);
    BOOST_CHECK_EQUAL(f0.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f0(0, 0), 1);
    BOOST_CHECK_EQUAL(f0(0, 1), 0);
    BOOST_CHECK_EQUAL(f0(1, 0), 3);
    BOOST_CHECK_EQUAL(f0(1, 1), 1);
    BOOST_CHECK_EQUAL(f0(2, 0), 5);
    BOOST_CHECK_EQUAL(f0(2, 1), 2);

    f0.get_underlying_policy().set_mapping(2, 2);
    BOOST_CHECK_EQUAL(f0.num_entries(), 6);
    BOOST_CHECK_EQUAL(f0.num_generators(), 3);
    BOOST_CHECK_EQUAL(f0.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f0(0, 0), 1);
    BOOST_CHECK_EQUAL(f0(0, 1), 2);
    BOOST_CHECK_EQUAL(f0(1, 0), 3);
    BOOST_CHECK_EQUAL(f0(1, 1), 4);
    BOOST_CHECK_EQUAL(f0(2, 0), 5);
    BOOST_CHECK_EQUAL(f0(2, 1), 6);

    f0.get_underlying_policy().set_mapping(9, -2);
    BOOST_CHECK_EQUAL(f0.num_entries(), 6);
    BOOST_CHECK_EQUAL(f0.num_generators(), 3);
    BOOST_CHECK_EQUAL(f0.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f0(0, 0), 1);
    BOOST_CHECK_EQUAL(f0(0, 1), 9);
    BOOST_CHECK_EQUAL(f0(1, 0), 3);
    BOOST_CHECK_EQUAL(f0(1, 1), 7);
    BOOST_CHECK_EQUAL(f0(2, 0), 5);
    BOOST_CHECK_EQUAL(f0(2, 1), 5);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_constructors, T, list_of_tested_variants) {
  test_constructors<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>();
  test_constructors<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>();
}

template <class F, typename T>
void test_utilities(T shift, T step) {
  F f0({1, 0, 2});
  f0.get_underlying_policy().set_mapping(shift, step);
  bool test = std::is_same_v<decltype(f0(0, 0)), T&>;
  BOOST_CHECK(test);

  BOOST_CHECK(!f0.is_plus_inf());
  BOOST_CHECK(!f0.is_minus_inf());
  BOOST_CHECK(!f0.is_nan());
  BOOST_CHECK(f0.is_finite());

  F f3;
  f3.get_underlying_policy().set_mapping(shift, step);
  BOOST_CHECK(!f3.is_plus_inf());
  BOOST_CHECK(f3.is_minus_inf());
  BOOST_CHECK(!f3.is_nan());
  BOOST_CHECK(!f3.is_finite());

  F f4 = F::minus_inf(2);
  f4.get_underlying_policy().set_mapping(shift, step);
  BOOST_CHECK(!f3.is_plus_inf());
  BOOST_CHECK(f3.is_minus_inf());
  BOOST_CHECK(!f3.is_nan());
  BOOST_CHECK(!f3.is_finite());

  F f5 = F::inf(2);
  f5.get_underlying_policy().set_mapping(shift, step);
  BOOST_CHECK(f5.is_plus_inf());
  BOOST_CHECK(!f5.is_minus_inf());
  BOOST_CHECK(!f5.is_nan());
  BOOST_CHECK(!f5.is_finite());

  if constexpr (std::numeric_limits<F>::has_quiet_NaN) {
    F f6 = F::nan(2);
    f6.get_underlying_policy().set_mapping(shift, step);
    BOOST_CHECK(!f6.is_plus_inf());
    BOOST_CHECK(!f6.is_minus_inf());
    BOOST_CHECK(f6.is_nan());
    BOOST_CHECK(!f6.is_finite());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_utilities, T, list_of_tested_variants) {
  test_utilities<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(0, 1);
  test_utilities<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(0, 1);

  test_utilities<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(-2, 2);
  test_utilities<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(0, 1);

  test_utilities<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(9, -2);
  test_utilities<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(0, 1);
}

template <class F, typename T>
void test_comparators() {
  const int num_param = 2;
  std::vector<T> v1, v2, v3, v4;

  if constexpr (F::ensures_1_criticality()) {
    v1 = {1, 0};
    v2 = {0, 0};
    v3 = {2, 0};
    v4 = {1, 0};
  } else {
    v1 = {4, 0, 3, 1, 2, 2};
    v2 = {0, 0, -2, 1, -1, 2, 0, 3};
    v3 = {5, 0};
    v4 = {5, 0, 3, 1, 1, 2};
  }

  F f1(v1.begin(), v1.end(), num_param);
  F f2(v2.begin(), v2.end(), num_param);
  F f3(v3.begin(), v3.end(), num_param);
  F f4(v4.begin(), v4.end(), num_param);
  F f5 = f1;
  f5.get_underlying_policy().set_mapping(0, 2);
  F f6 = f1;
  f6.get_underlying_policy().set_mapping(-9, 2);

  BOOST_CHECK(!(f1 < f1));
  BOOST_CHECK(!(f1 < f2));
  BOOST_CHECK(f1 < f3);
  BOOST_CHECK(!(f1 < f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 < F::nan(num_param)));
  BOOST_CHECK(f1 < F::inf(num_param));
  BOOST_CHECK(!(f1 < F::minus_inf(num_param)));
  BOOST_CHECK(!(f1 < f5));
  BOOST_CHECK(!(f1 < f6));

  BOOST_CHECK(f1 <= f1);
  BOOST_CHECK(!(f1 <= f2));
  BOOST_CHECK(f1 <= f3);
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(f1 <= f4);
  else
    BOOST_CHECK(!(f1 <= f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 <= F::nan(num_param)));
  BOOST_CHECK(f1 <= F::inf(num_param));
  BOOST_CHECK(!(f1 <= F::minus_inf(num_param)));
  BOOST_CHECK(f1 <= f5);
  BOOST_CHECK(!(f1 <= f6));

  BOOST_CHECK(!(f1 > f1));
  BOOST_CHECK(f1 > f2);
  BOOST_CHECK(!(f1 > f3));
  BOOST_CHECK(!(f1 > f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 > F::nan(num_param)));
  BOOST_CHECK(!(f1 > F::inf(num_param)));
  BOOST_CHECK(f1 > F::minus_inf(num_param));
  BOOST_CHECK(!(f1 > f5));
  BOOST_CHECK(f1 > f6);

  BOOST_CHECK(f1 >= f1);
  BOOST_CHECK(f1 >= f2);
  BOOST_CHECK(!(f1 >= f3));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(f1 >= f4);
  else
    BOOST_CHECK(!(f1 >= f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 >= F::nan(num_param)));
  BOOST_CHECK(!(f1 >= F::inf(num_param)));
  BOOST_CHECK(f1 >= F::minus_inf(num_param));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(f1 >= f5);
  else
    BOOST_CHECK(!(f1 >= f5));
  BOOST_CHECK(f1 >= f6);

  BOOST_CHECK(f1 == f1);
  BOOST_CHECK(!(f1 == f2));
  BOOST_CHECK(!(f1 == f3));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(f1 == f4);
  else
    BOOST_CHECK(!(f1 == f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(!(f1 == F::nan(num_param)));
  BOOST_CHECK(!(f1 == F::inf(num_param)));
  BOOST_CHECK(!(f1 == F::minus_inf(num_param)));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(f1 == f5);
  else
    BOOST_CHECK(!(f1 == f5));
  BOOST_CHECK(!(f1 == f6));

  BOOST_CHECK(!(f1 != f1));
  BOOST_CHECK(f1 != f2);
  BOOST_CHECK(f1 != f3);
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(!(f1 != f4));
  else
    BOOST_CHECK(f1 != f4);
  if constexpr (std::numeric_limits<F>::has_quiet_NaN) BOOST_CHECK(f1 != F::nan(num_param));
  BOOST_CHECK(f1 != F::inf(num_param));
  BOOST_CHECK(f1 != F::minus_inf(num_param));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(!(f1 != f5));
  else
    BOOST_CHECK(f1 != f5);
  BOOST_CHECK(f1 != f6);
}

template <class F, typename T>
void test_lex_comparators() {
  const int num_param = 2;
  std::vector<T> v1, v2, v3, v4, v5;

  if constexpr (F::ensures_1_criticality()) {
    v2 = {0, 0};
    v3 = {1, 0};
    v1 = {2, 0};
    v4 = {3, 0};
    v5 = {3, 0};
  } else {
    v2 = {0, 0, -2, 1, -1, 2, 0, 3};
    v3 = {2, 0};
    v1 = {2, 0, 3, 1, 2, 2};
    v4 = {3, 0, 3, 1, 1, 2};
    v5 = {3, 0, 3, 1, 1, 2};
  }

  F f1(v1.begin(), v1.end(), num_param);
  F f2(v2.begin(), v2.end(), num_param);
  F f3(v3.begin(), v3.end(), num_param);
  F f4(v4.begin(), v4.end(), num_param);
  F f5(v5.begin(), v5.end(), num_param);
  F f6 = f1;
  f6.get_underlying_policy().set_mapping(0, 2);
  F f7 = f1;
  f7.get_underlying_policy().set_mapping(-9, 2);
  F f8 = f1;
  f8.get_underlying_policy().set_mapping(3, -2);

  BOOST_CHECK(is_less_or_equal_than_lexicographically(f1, f1));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f1, f2));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f1, f3));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f1, f4));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f1, f6));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f1, f7));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f1, f8));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_less_or_equal_than_lexicographically(f1, F::nan(num_param)));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f1, F::inf(num_param)));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f1, F::minus_inf(num_param)));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(is_less_or_equal_than_lexicographically(f6, f1));
  else
    BOOST_CHECK(!is_less_or_equal_than_lexicographically(f6, f1));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f7, f1));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f8, f1));

  BOOST_CHECK(!is_strict_less_than_lexicographically(f1, f1));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f1, f2));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f1, f3));
  BOOST_CHECK(is_strict_less_than_lexicographically(f1, f4));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(!is_strict_less_than_lexicographically(f1, f6));
  else
    BOOST_CHECK(is_strict_less_than_lexicographically(f1, f6));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f1, f7));
  BOOST_CHECK(is_strict_less_than_lexicographically(f1, f8));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_strict_less_than_lexicographically(f1, F::nan(num_param)));
  BOOST_CHECK(is_strict_less_than_lexicographically(f1, F::inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f1, F::minus_inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f6, f1));
  BOOST_CHECK(is_strict_less_than_lexicographically(f7, f1));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f8, f1));

  BOOST_CHECK(is_less_or_equal_than_lexicographically(f2, f1));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f2, f2));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f2, f3));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f2, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_less_or_equal_than_lexicographically(f2, F::nan(num_param)));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f2, F::inf(num_param)));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f2, F::minus_inf(num_param)));

  BOOST_CHECK(is_strict_less_than_lexicographically(f2, f1));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f2, f2));
  BOOST_CHECK(is_strict_less_than_lexicographically(f2, f3));
  BOOST_CHECK(is_strict_less_than_lexicographically(f2, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_strict_less_than_lexicographically(f2, F::nan(num_param)));
  BOOST_CHECK(is_strict_less_than_lexicographically(f2, F::inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f2, F::minus_inf(num_param)));

  BOOST_CHECK(is_less_or_equal_than_lexicographically(f3, f1));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f3, f2));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f3, f3));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f3, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_less_or_equal_than_lexicographically(f3, F::nan(num_param)));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f3, F::inf(num_param)));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f3, F::minus_inf(num_param)));

  BOOST_CHECK(is_strict_less_than_lexicographically(f3, f1));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f3, f2));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f3, f3));
  BOOST_CHECK(is_strict_less_than_lexicographically(f3, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_strict_less_than_lexicographically(f3, F::nan(num_param)));
  BOOST_CHECK(is_strict_less_than_lexicographically(f3, F::inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f3, F::minus_inf(num_param)));

  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f4, f1));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f4, f2));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f4, f3));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f4, f4));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f4, f5));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f5, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_less_or_equal_than_lexicographically(f4, F::nan(num_param)));
  BOOST_CHECK(is_less_or_equal_than_lexicographically(f4, F::inf(num_param)));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically(f4, F::minus_inf(num_param)));

  BOOST_CHECK(!is_strict_less_than_lexicographically(f4, f1));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f4, f2));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f4, f3));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f4, f4));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f4, f5));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f5, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_strict_less_than_lexicographically(f4, F::nan(num_param)));
  BOOST_CHECK(is_strict_less_than_lexicographically(f4, F::inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically(f4, F::minus_inf(num_param)));
}

template <class F, typename T>
void test_co_lex_comparators() {
  using namespace Gudhi::multi_filtration;

  const int num_param = 2;
  std::vector<T> v1, v2, v3, v4, v5;

  if constexpr (F::ensures_1_criticality()) {
    v1 = {2, 0};
    v2 = {4, 0};
    v3 = {1, 0};
    v4 = {3, 0};
    v5 = {3, 0};
  } else {
    v1 = {2, 0, 3, 1, 2, 2};
    v2 = {0, 0, -2, 1, -1, 2, 0, 3};
    v3 = {2, 0};
    v4 = {3, 0, 3, 1, 3, 2};
    v5 = {3, 0, 3, 1, 3, 2};
  }

  F f1(v1.begin(), v1.end(), num_param);
  F f2(v2.begin(), v2.end(), num_param);
  F f3(v3.begin(), v3.end(), num_param);
  F f4(v4.begin(), v4.end(), num_param);
  F f5(v5.begin(), v5.end(), num_param);
  F f6 = f1;
  f6.get_underlying_policy().set_mapping(0, 2);
  F f7 = f1;
  f7.get_underlying_policy().set_mapping(-9, 2);
  F f8 = f1;
  f8.get_underlying_policy().set_mapping(3, -2);

  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f1, f1));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f1, f2));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f1, f3));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f1, f4));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f1, f6));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f1, f7));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f1, f8));
  else
    BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f1, f8));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f1, F::nan(num_param)));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f1, F::inf(num_param)));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f1, F::minus_inf(num_param)));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f6, f1));
  else
    BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f6, f1));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f7, f1));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f8, f1));
  else
    BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f8, f1));

  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f1, f1));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f1, f2));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f1, f3));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f1, f4));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f1, f6));
  else
    BOOST_CHECK(is_strict_less_than_lexicographically<true>(f1, f6));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f1, f7));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f1, f8));
  else
    BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f1, f8));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_strict_less_than_lexicographically<true>(f1, F::nan(num_param)));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f1, F::inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f1, F::minus_inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f6, f1));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f7, f1));
  if constexpr (F::ensures_1_criticality())
    BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f8, f1));
  else
    BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f8, f1));

  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f2, f1));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f2, f2));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f2, f3));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f2, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f2, F::nan(num_param)));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f2, F::inf(num_param)));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f2, F::minus_inf(num_param)));

  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f2, f1));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f2, f2));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f2, f3));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f2, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_strict_less_than_lexicographically<true>(f2, F::nan(num_param)));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f2, F::inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f2, F::minus_inf(num_param)));

  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f3, f1));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f3, f2));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f3, f3));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f3, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f3, F::nan(num_param)));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f3, F::inf(num_param)));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f3, F::minus_inf(num_param)));

  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f3, f1));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f3, f2));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f3, f3));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f3, f4));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_strict_less_than_lexicographically<true>(f3, F::nan(num_param)));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f3, F::inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f3, F::minus_inf(num_param)));

  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f4, f1));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f4, f2));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f4, f3));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f4, f4));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f5, f4));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f4, f5));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f4, F::nan(num_param)));
  BOOST_CHECK(is_less_or_equal_than_lexicographically<true>(f4, F::inf(num_param)));
  BOOST_CHECK(!is_less_or_equal_than_lexicographically<true>(f4, F::minus_inf(num_param)));

  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f4, f1));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f4, f2));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f4, f3));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f4, f4));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f5, f4));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f4, f5));
  if constexpr (std::numeric_limits<F>::has_quiet_NaN)
    BOOST_CHECK(is_strict_less_than_lexicographically<true>(f4, F::nan(num_param)));
  BOOST_CHECK(is_strict_less_than_lexicographically<true>(f4, F::inf(num_param)));
  BOOST_CHECK(!is_strict_less_than_lexicographically<true>(f4, F::minus_inf(num_param)));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_comparators, T, list_of_tested_variants) {
  test_comparators<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>();
  test_comparators<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>();
  test_lex_comparators<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>();
  test_lex_comparators<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>();
  test_co_lex_comparators<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>();
  test_co_lex_comparators<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>();
}

template <class F, typename T>
void test_operators() {
  const int num_param = 2;

  F f({-1, 0});
  F f2({2, 0});
  F f3({F::T_m_inf, 0});
  F f4({F::T_inf, 0});
  // TODO: tests with more than 1 generator

  F res = -f;
  BOOST_CHECK_EQUAL(res(0, 0), 1);
  BOOST_CHECK_EQUAL(res(0, 1), 0);
  BOOST_CHECK((-F::inf(num_param)).is_minus_inf());
  BOOST_CHECK((-F::minus_inf(num_param)).is_plus_inf());
  BOOST_CHECK((-F::nan(num_param)).is_nan());

  res = f - f2;
  BOOST_CHECK_EQUAL(res(0, 0), -3);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f - f3;
  BOOST_CHECK_EQUAL(res(0, 0), f4(0, 0));
  BOOST_CHECK_EQUAL(res(0, 1), f4(0, 1));

  res = f3 - f;
  BOOST_CHECK_EQUAL(res(0, 0), f3(0, 0));
  BOOST_CHECK_EQUAL(res(0, 1), f3(0, 1));

  res = T(5) - f;
  BOOST_CHECK_EQUAL(res(0, 0), 6);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f - T(5);
  BOOST_CHECK_EQUAL(res(0, 0), -6);
  BOOST_CHECK_EQUAL(res(0, 1), 0);
  BOOST_CHECK((f - F::inf(num_param)).is_minus_inf());
  BOOST_CHECK((F::inf(num_param) - f).is_plus_inf());
  BOOST_CHECK((f - F::minus_inf(num_param)).is_plus_inf());
  BOOST_CHECK((F::minus_inf(num_param) - f).is_minus_inf());
  BOOST_CHECK((f - F::nan(num_param)).is_nan());
  BOOST_CHECK((F::nan(num_param) - f).is_nan());

  res = f3 - f3;
  BOOST_CHECK(res.is_nan());

  res = f3 - f4;
  BOOST_CHECK_EQUAL(res(0, 0), f3(0, 0));
  BOOST_CHECK_EQUAL(res(0, 1), f3(0, 1));

  res = f + f2;
  BOOST_CHECK_EQUAL(res(0, 0), 1);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f + f3;
  BOOST_CHECK_EQUAL(res(0, 0), f3(0, 0));
  BOOST_CHECK_EQUAL(res(0, 1), f3(0, 1));

  res = T(5) + f;
  BOOST_CHECK_EQUAL(res(0, 0), 4);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f + T(5);
  BOOST_CHECK_EQUAL(res(0, 0), 4);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  BOOST_CHECK((f + F::inf(num_param)).is_plus_inf());
  BOOST_CHECK((F::inf(num_param) + f).is_plus_inf());
  BOOST_CHECK((f + F::minus_inf(num_param)).is_minus_inf());
  BOOST_CHECK((F::minus_inf(num_param) + f).is_minus_inf());
  BOOST_CHECK((f + F::nan(num_param)).is_nan());
  BOOST_CHECK((F::nan(num_param) + f).is_nan());

  res = f3 + f4;
  BOOST_CHECK(res.is_nan());

  res = f3 + f3;
  BOOST_CHECK_EQUAL(res(0, 0), f3(0, 0));
  BOOST_CHECK_EQUAL(res(0, 1), f3(0, 1));

  res = f * f2;
  BOOST_CHECK_EQUAL(res(0, 0), -2);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f * f3;
  BOOST_CHECK_EQUAL(res(0, 0), f4(0, 0));
  BOOST_CHECK_EQUAL(res(0, 1), f4(0, 1));

  res = T(0) * f3;
  BOOST_CHECK(res.is_nan());

  res = T(5) * f;
  BOOST_CHECK_EQUAL(res(0, 0), -5);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f * T(5);
  BOOST_CHECK_EQUAL(res(0, 0), -5);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f * F::inf(num_param);
  BOOST_CHECK_EQUAL(res(0, 0), F::T_m_inf);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = F::inf(num_param) * f;
  BOOST_CHECK_EQUAL(res(0, 0), F::T_m_inf);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f * F::minus_inf(num_param);
  BOOST_CHECK_EQUAL(res(0, 0), F::T_inf);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = F::minus_inf(num_param) * f;
  BOOST_CHECK_EQUAL(res(0, 0), F::T_inf);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  BOOST_CHECK((f * F::nan(num_param)).is_nan());
  BOOST_CHECK((F::nan(num_param) * f).is_nan());

  res = f3 * f3;
  BOOST_CHECK(res.is_plus_inf());
  res = f3 * f4;
  BOOST_CHECK(res.is_minus_inf());

  res = f / f2;
  BOOST_CHECK_EQUAL(res(0, 0), static_cast<T>(-0.5));
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f / f3;
  BOOST_CHECK_EQUAL(res(0, 0), 0);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f3 / f;
  BOOST_CHECK_EQUAL(res(0, 0), f4(0, 0));
  BOOST_CHECK_EQUAL(res(0, 1), f4(0, 1));

  res = f3 / T(0);
  BOOST_CHECK(res.is_nan());

  res = T(5) / f;
  BOOST_CHECK_EQUAL(res(0, 0), -5);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f / T(5);
  BOOST_CHECK_EQUAL(res(0, 0), static_cast<T>(-0.2));
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = f / F::inf(num_param);
  BOOST_CHECK_EQUAL(res(0, 0), 0);
  BOOST_CHECK_EQUAL(res(0, 1), 0);
  res = F::inf(num_param) / f;
  BOOST_CHECK_EQUAL(res(0, 0), f3(0, 0));
  BOOST_CHECK_EQUAL(res(0, 1), f3(0, 1));

  res = f / F::minus_inf(num_param);
  BOOST_CHECK_EQUAL(res(0, 0), 0);
  BOOST_CHECK_EQUAL(res(0, 1), 0);
  res = F::minus_inf(num_param) / f;
  BOOST_CHECK_EQUAL(res(0, 0), f4(0, 0));
  BOOST_CHECK_EQUAL(res(0, 1), f4(0, 1));

  res = f / F::nan(num_param);
  BOOST_CHECK(res.is_nan());
  res = F::nan(num_param) / f;
  BOOST_CHECK(res.is_nan());

  res = f3 / f3;
  BOOST_CHECK(res.is_nan());
  res = f3 / f4;
  BOOST_CHECK(res.is_nan());
  res = f / F({0, 0, 0});
  BOOST_CHECK(res.is_nan());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_operators, T, list_of_tested_variants) {
  test_operators<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>();
  test_operators<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>();
}

template <class F, typename T>
void test_modifiers1(T shift, T step) {
  const int num_param = 2;
  const T inf = F::T_inf;
  const T m_inf = F::T_m_inf;
  std::vector<T> v;

  auto get_value = [shift, step](T v) -> T { return shift + (v * step); };

  if constexpr (F::ensures_1_criticality()) {
    v = {5, 0};
  } else {
    v = {5, 0, 3, 1, 1, 2, 7, 3};
  }

  F f1(v.begin(), v.end(), num_param);
  F f2(v.begin(), v.end(), num_param);
  F f3(v.begin(), v.end(), num_param);
  F f4(v.begin(), v.end(), num_param);
  F f5(v.begin(), v.end(), num_param);
  F f6(v.begin(), v.end(), num_param);
  F f7(v.begin(), v.end(), num_param);
  F f8(v.begin(), v.end(), num_param);
  f1.get_underlying_policy().set_mapping(shift, step);
  f2.get_underlying_policy().set_mapping(shift, step);
  f3.get_underlying_policy().set_mapping(shift, step);
  f4.get_underlying_policy().set_mapping(shift, step);
  f5.get_underlying_policy().set_mapping(shift, step);
  f6.get_underlying_policy().set_mapping(shift, step);
  f7.get_underlying_policy().set_mapping(shift, step);
  f8.get_underlying_policy().set_mapping(shift, step);

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f1(0, 0), 5);
    BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
  } else {
    BOOST_CHECK_EQUAL(f1(0, 0), 5);
    BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f1(1, 0), 3);
    BOOST_CHECK_EQUAL(f1(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f1(2, 0), 1);
    BOOST_CHECK_EQUAL(f1(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f1(3, 0), 7);
    BOOST_CHECK_EQUAL(f1(3, 1), get_value(3));
  }

  f1.push_to_least_common_upper_bound(F::minus_inf(num_param));
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f1(0, 0), 5);
    BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
  } else {
    BOOST_CHECK_EQUAL(f1(0, 0), 5);
    BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f1(1, 0), 3);
    BOOST_CHECK_EQUAL(f1(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f1(2, 0), 1);
    BOOST_CHECK_EQUAL(f1(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f1(3, 0), 7);
    BOOST_CHECK_EQUAL(f1(3, 1), get_value(3));
  }

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_THROW((f2.push_to_least_common_upper_bound({0, get_value(1)})), std::invalid_argument);
  } else {
    f2.push_to_least_common_upper_bound({0, get_value(1)});
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f2(0, 0), m_inf);
      BOOST_CHECK_EQUAL(f2(1, 0), 5);
    } else {
      BOOST_CHECK_EQUAL(f2(0, 0), inf);
      BOOST_CHECK_EQUAL(f2(1, 0), 3);
    }
    BOOST_CHECK_EQUAL(f2(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f2(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f2(2, 0), 1);
    BOOST_CHECK_EQUAL(f2(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f2(3, 0), 7);
    BOOST_CHECK_EQUAL(f2(3, 1), get_value(3));
  }

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_THROW((f3.push_to_least_common_upper_bound({0, get_value(4)})), std::invalid_argument);
  } else {
    f3.push_to_least_common_upper_bound({0, get_value(4)});
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f3(0, 0), m_inf);
      BOOST_CHECK_EQUAL(f3(1, 0), m_inf);
      BOOST_CHECK_EQUAL(f3(2, 0), m_inf);
      BOOST_CHECK_EQUAL(f3(3, 0), m_inf);
      BOOST_CHECK_EQUAL(f3(4, 0), 7);
    } else {
      BOOST_CHECK_EQUAL(f3(0, 0), inf);
      BOOST_CHECK_EQUAL(f3(1, 0), inf);
      BOOST_CHECK_EQUAL(f3(2, 0), inf);
      BOOST_CHECK_EQUAL(f3(3, 0), inf);
      BOOST_CHECK_EQUAL(f3(4, 0), 1);
    }
    BOOST_CHECK_EQUAL(f3(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f3(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f3(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f3(3, 1), get_value(3));
    BOOST_CHECK_EQUAL(f3(4, 1), get_value(4));
  }

  if constexpr (F::ensures_1_criticality()) {
    f4.push_to_least_common_upper_bound({2, get_value(0)});
    BOOST_CHECK_EQUAL(f4(0, 0), 5);
    BOOST_CHECK_EQUAL(f4(0, 1), get_value(0));
  } else {
    f4.push_to_least_common_upper_bound({2, get_value(1)});
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f4(0, 0), m_inf);
      BOOST_CHECK_EQUAL(f4(1, 0), 5);
    } else {
      BOOST_CHECK_EQUAL(f4(0, 0), inf);
      BOOST_CHECK_EQUAL(f4(1, 0), 3);
    }
    BOOST_CHECK_EQUAL(f4(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f4(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f4(2, 0), 2);
    BOOST_CHECK_EQUAL(f4(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f4(3, 0), 7);
    BOOST_CHECK_EQUAL(f4(3, 1), get_value(3));
  }

  if constexpr (F::ensures_1_criticality()) {
    f5.push_to_least_common_upper_bound({4, get_value(0)});
    BOOST_CHECK_EQUAL(f5(0, 0), 5);
    BOOST_CHECK_EQUAL(f5(0, 1), get_value(0));
  } else {
    f5.push_to_least_common_upper_bound({4, get_value(1)});
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f5(0, 0), m_inf);
      BOOST_CHECK_EQUAL(f5(1, 0), 5);
    } else {
      BOOST_CHECK_EQUAL(f5(0, 0), inf);
      BOOST_CHECK_EQUAL(f5(1, 0), 4);
    }
    BOOST_CHECK_EQUAL(f5(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f5(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f5(2, 0), 4);
    BOOST_CHECK_EQUAL(f5(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f5(3, 0), 7);
    BOOST_CHECK_EQUAL(f5(3, 1), get_value(3));
  }

  f6.push_to_least_common_upper_bound({inf, get_value(0)});
  if constexpr (F::ensures_1_criticality()) {
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f6(0, 0), inf);
      BOOST_CHECK_EQUAL(f6(0, 1), get_value(0));
    } else {
      BOOST_CHECK(f6.is_plus_inf());
    }
  } else {
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f6(0, 0), inf);
      BOOST_CHECK_EQUAL(f6(0, 1), get_value(0));
      BOOST_CHECK_EQUAL(f6(1, 0), inf);
      BOOST_CHECK_EQUAL(f6(1, 1), get_value(1));
      BOOST_CHECK_EQUAL(f6(2, 0), inf);
      BOOST_CHECK_EQUAL(f6(2, 1), get_value(2));
      BOOST_CHECK_EQUAL(f6(3, 0), inf);
      BOOST_CHECK_EQUAL(f6(3, 1), get_value(3));
    } else {
      BOOST_CHECK(f6.is_plus_inf());
    }
  }

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_THROW((f7.push_to_least_common_upper_bound({9, get_value(2)})), std::invalid_argument);
  } else {
    f7.push_to_least_common_upper_bound({9, get_value(2)});
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f7(0, 0), m_inf);
      BOOST_CHECK_EQUAL(f7(1, 0), m_inf);
    } else {
      BOOST_CHECK_EQUAL(f7(0, 0), inf);
      BOOST_CHECK_EQUAL(f7(1, 0), inf);
    }
    BOOST_CHECK_EQUAL(f7(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f7(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f7(2, 0), 9);
    BOOST_CHECK_EQUAL(f7(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f7(3, 0), 9);
    BOOST_CHECK_EQUAL(f7(3, 1), get_value(3));
  }

  if constexpr (F::ensures_1_criticality()) {
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_THROW((f8.push_to_least_common_upper_bound({9, inf})), std::invalid_argument);
    } else {
      f8.push_to_least_common_upper_bound({9, inf});
      BOOST_CHECK(f8.is_plus_inf());
    }
  } else {
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_THROW((f8.push_to_least_common_upper_bound({9, inf})), std::invalid_argument);
    } else {
      f8.push_to_least_common_upper_bound({9, inf});
      BOOST_CHECK(f8.is_plus_inf());
    }
  }

  f1.push_to_least_common_upper_bound(F::nan(num_param));
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f1(0, 0), 5);
    BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
  } else {
    BOOST_CHECK_EQUAL(f1(0, 0), 5);
    BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f1(1, 0), 3);
    BOOST_CHECK_EQUAL(f1(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f1(2, 0), 1);
    BOOST_CHECK_EQUAL(f1(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f1(3, 0), 7);
    BOOST_CHECK_EQUAL(f1(3, 1), get_value(3));
  }

  F b = F::minus_inf(num_param);
  b.get_underlying_policy().set_mapping(shift, step);
  b.push_to_least_common_upper_bound({3, get_value(0)});
  BOOST_CHECK_EQUAL(b(0, 0), 3);
  BOOST_CHECK_EQUAL(b(0, 1), get_value(0));
}

template <class F, typename T>
void test_modifiers2(T shift, T step) {
  const int num_param = 2;
  const T inf = F::T_inf;
  std::vector<T> v;

  auto get_value = [shift, step](T v) -> T { return shift + (v * step); };

  if constexpr (F::ensures_1_criticality()) {
    v = {5, 0};
  } else {
    v = {5, 0, 3, 1, 1, 2, 7, 3};
  }

  F f1(v.begin(), v.end(), num_param);
  F f2(v.begin(), v.end(), num_param);
  F f3(v.begin(), v.end(), num_param);
  F f4(v.begin(), v.end(), num_param);
  F f5(v.begin(), v.end(), num_param);
  F f6(v.begin(), v.end(), num_param);
  F f7(v.begin(), v.end(), num_param);
  F f8(v.begin(), v.end(), num_param);
  f1.get_underlying_policy().set_mapping(shift, step);
  f2.get_underlying_policy().set_mapping(shift, step);
  f3.get_underlying_policy().set_mapping(shift, step);
  f4.get_underlying_policy().set_mapping(shift, step);
  f5.get_underlying_policy().set_mapping(shift, step);
  f6.get_underlying_policy().set_mapping(shift, step);
  f7.get_underlying_policy().set_mapping(shift, step);
  f8.get_underlying_policy().set_mapping(shift, step);

  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f1(0, 0), 5);
    BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
  } else {
    BOOST_CHECK_EQUAL(f1(0, 0), 5);
    BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f1(1, 0), 3);
    BOOST_CHECK_EQUAL(f1(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f1(2, 0), 1);
    BOOST_CHECK_EQUAL(f1(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f1(3, 0), 7);
    BOOST_CHECK_EQUAL(f1(3, 1), get_value(3));
  }

  f1.pull_to_greatest_common_lower_bound(F::minus_inf(num_param));
  BOOST_CHECK(f1.is_minus_inf());

  f2.pull_to_greatest_common_lower_bound({0, get_value(1)});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f2(0, 0), 0);
    BOOST_CHECK_EQUAL(f2(0, 1), get_value(0));
  } else {
    BOOST_CHECK_EQUAL(f2(0, 0), 0);
    BOOST_CHECK_EQUAL(f2(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f2(1, 0), 0);
    BOOST_CHECK_EQUAL(f2(1, 1), get_value(1));
  }

  f3.pull_to_greatest_common_lower_bound({0, get_value(4)});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f3(0, 0), 0);
    BOOST_CHECK_EQUAL(f3(0, 1), get_value(0));
  } else {
    BOOST_CHECK_EQUAL(f3(0, 0), 0);
    BOOST_CHECK_EQUAL(f3(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f3(1, 0), 0);
    BOOST_CHECK_EQUAL(f3(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f3(2, 0), 0);
    BOOST_CHECK_EQUAL(f3(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f3(3, 0), 0);
    BOOST_CHECK_EQUAL(f3(3, 1), get_value(3));
  }

  f4.pull_to_greatest_common_lower_bound({2, get_value(1)});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f4(0, 0), 2);
    BOOST_CHECK_EQUAL(f4(0, 1), get_value(0));
  } else {
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f4(1, 0), 2);
    } else {
      BOOST_CHECK_EQUAL(f4(1, 0), 1);
    }
    BOOST_CHECK_EQUAL(f4(0, 0), 2);
    BOOST_CHECK_EQUAL(f4(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f4(1, 1), get_value(1));
  }

  f5.pull_to_greatest_common_lower_bound({4, get_value(1)});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f5(0, 0), 4);
    BOOST_CHECK_EQUAL(f5(0, 1), get_value(0));
  } else {
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f5(1, 0), 4);
    } else {
      BOOST_CHECK_EQUAL(f5(1, 0), 1);
    }
    BOOST_CHECK_EQUAL(f5(0, 0), 4);
    BOOST_CHECK_EQUAL(f5(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f5(1, 1), get_value(1));
  }

  f6.pull_to_greatest_common_lower_bound({inf, get_value(0)});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f6(0, 0), 5);
    BOOST_CHECK_EQUAL(f6(0, 1), get_value(0));
  } else {
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f6(0, 0), 7);
    } else {
      BOOST_CHECK_EQUAL(f6(0, 0), 1);
    }
    BOOST_CHECK_EQUAL(f6(0, 1), get_value(0));
  }

  f7.pull_to_greatest_common_lower_bound({9, get_value(2)});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f7(0, 0), 5);
    BOOST_CHECK_EQUAL(f7(0, 1), get_value(0));
  } else {
    if constexpr (F::has_negative_cones()) {
      BOOST_CHECK_EQUAL(f7(2, 0), 7);
    } else {
      BOOST_CHECK_EQUAL(f7(2, 0), 1);
    }
    BOOST_CHECK_EQUAL(f7(0, 0), 5);
    BOOST_CHECK_EQUAL(f7(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f7(1, 0), 3);
    BOOST_CHECK_EQUAL(f7(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f7(2, 1), get_value(2));
  }

  f8.pull_to_greatest_common_lower_bound({9, inf});
  if constexpr (F::ensures_1_criticality()) {
    BOOST_CHECK_EQUAL(f8(0, 0), 5);
    BOOST_CHECK_EQUAL(f8(0, 1), get_value(0));
  } else {
    BOOST_CHECK_EQUAL(f8(0, 0), 5);
    BOOST_CHECK_EQUAL(f8(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f8(1, 0), 3);
    BOOST_CHECK_EQUAL(f8(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f8(2, 0), 1);
    BOOST_CHECK_EQUAL(f8(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f8(3, 0), 7);
    BOOST_CHECK_EQUAL(f8(3, 1), get_value(3));
  }

  f1.pull_to_greatest_common_lower_bound(F::nan(num_param));
  BOOST_CHECK(f1.is_minus_inf());

  if constexpr (!F::has_negative_cones()) {
    F a = F::inf(num_param);
    a.get_underlying_policy().set_mapping(shift, step);
    a.pull_to_greatest_common_lower_bound({3, get_value(0)});
    BOOST_CHECK_EQUAL(a(0, 0), 3);
    BOOST_CHECK_EQUAL(a(0, 1), get_value(0));
  }
}

template <class F, typename T>
void test_modifiers3(T shift, T step) {
  F f1({7, 0});
  f1.get_underlying_policy().set_mapping(shift, step);

  auto get_value = [shift, step](T v) -> T { return shift + (v * step); };

  std::vector<std::vector<int>> grid = {{0, 3, 6, 9}, {0, 1, 2, 3}, {0, 4, 8, 16}};

  f1.project_onto_grid(grid, true);
  BOOST_CHECK_EQUAL(f1(0, 0), 2);
  BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));

  f1 = F({7, 0});
  f1.get_underlying_policy().set_mapping(shift, step);
  f1.project_onto_grid(grid, false);
  BOOST_CHECK_EQUAL(f1(0, 0), 6);
  BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));

  if constexpr (!F::ensures_1_criticality()) {
    T def = F::has_negative_cones() ? F::T_m_inf : F::T_inf;
    f1.set_num_generators(5);
    BOOST_CHECK_EQUAL(f1.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f1.num_generators(), 5);
    BOOST_CHECK_EQUAL(f1.num_entries(), 10);
    BOOST_CHECK_EQUAL(f1(0, 0), 6);
    BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f1(1, 0), def);
    BOOST_CHECK_EQUAL(f1(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f1(2, 0), def);
    BOOST_CHECK_EQUAL(f1(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f1(3, 0), def);
    BOOST_CHECK_EQUAL(f1(3, 1), get_value(3));
    BOOST_CHECK_EQUAL(f1(4, 0), def);
    BOOST_CHECK_EQUAL(f1(4, 1), get_value(4));
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_modifiers, T, list_of_tested_variants) {
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(0, 1);
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(0, 1);
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>(0, 1);
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>(0, 1);

  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(-2, 2);
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(-2, 2);
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>(-2, 2);
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>(-2, 2);

  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(9, -2);
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(9, -2);
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>(9, -2);
  test_modifiers1<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>(9, -2);

  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(0, 1);
  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(0, 1);
  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>(0, 1);
  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>(0, 1);

  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(-2, 2);
  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(-2, 2);
  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>(-2, 2);
  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>(-2, 2);

  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(9, -2);
  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(9, -2);
  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>(9, -2);
  test_modifiers2<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>(9, -2);

  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(0, 1);
  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(0, 1);
  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>(0, 1);
  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>(0, 1);

  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(-2, 2);
  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(-2, 2);
  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>(-2, 2);
  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>(-2, 2);

  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(9, -2);
  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>(9, -2);
  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>(9, -2);
  test_modifiers3<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>(9, -2);
}

template <class F, typename T>
void test_add_generators(T shift, T step) {
  const int num_param = 2;

  auto get_value = [shift, step](T v) -> T { return shift + (v * step); };

  F f({1, 0});
  f.get_underlying_policy().set_mapping(shift, step);

  BOOST_CHECK_EQUAL(f.num_generators(), 1);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), 1);
  BOOST_CHECK_EQUAL(f(0, 1), get_value(0));

  bool res = f.add_generator({2, get_value(1)});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), 1);
  BOOST_CHECK_EQUAL(f(0, 1), get_value(0));
  BOOST_CHECK_EQUAL(f(1, 0), 2);
  BOOST_CHECK_EQUAL(f(1, 1), get_value(1));

  res = f.add_generator({-2, get_value(0)});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), -2);
  BOOST_CHECK_EQUAL(f(0, 1), get_value(0));
  BOOST_CHECK_EQUAL(f(1, 0), 2);
  BOOST_CHECK_EQUAL(f(1, 1), get_value(1));

  res = f.add_generator({3, get_value(0)});
  BOOST_CHECK(!res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), -2);
  BOOST_CHECK_EQUAL(f(0, 1), get_value(0));
  BOOST_CHECK_EQUAL(f(1, 0), 2);
  BOOST_CHECK_EQUAL(f(1, 1), get_value(1));

  res = f.add_generator({F::T_inf, get_value(0)});
  BOOST_CHECK(!res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), -2);
  BOOST_CHECK_EQUAL(f(0, 1), get_value(0));
  BOOST_CHECK_EQUAL(f(1, 0), 2);
  BOOST_CHECK_EQUAL(f(1, 1), get_value(1));

  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    res = f.add_generator({std::numeric_limits<T>::quiet_NaN(), get_value(0)});
    BOOST_CHECK(!res);
    BOOST_CHECK_EQUAL(f.num_generators(), 2);
    BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
    BOOST_CHECK_EQUAL(f(0, 0), -2);
    BOOST_CHECK_EQUAL(f(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f(1, 0), 2);
    BOOST_CHECK_EQUAL(f(1, 1), get_value(1));
  }

  res = f.add_generator({F::T_m_inf, get_value(0)});
  BOOST_CHECK(res);
  BOOST_CHECK_EQUAL(f.num_generators(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), num_param);
  BOOST_CHECK_EQUAL(f(0, 0), F::T_m_inf);
  BOOST_CHECK_EQUAL(f(0, 1), get_value(0));
  BOOST_CHECK_EQUAL(f(1, 0), 2);
  BOOST_CHECK_EQUAL(f(1, 1), get_value(1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_add_generators, T, list_of_tested_variants) {
  test_add_generators<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(0, 1);
  test_add_generators<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(-2, 2);
  test_add_generators<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(9, -2);

  Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true> f({1, 0});
  BOOST_CHECK_THROW(f.add_generator({1, 1}), std::logic_error);
}

template <class F, typename T>
void test_friends() {
  T inf = F::has_negative_cones() ? F::T_m_inf : F::T_inf;
  F f({1, 2}, 2);

  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(6))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, F({4, 5, 3}, 2)), static_cast<T>(std::sqrt(T(2))));
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {3, 2, 5, 9}), 3);
  BOOST_CHECK_EQUAL(compute_norm<double>(f), std::sqrt(double(6.)));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to<double>(f, F({4, 5, 3}, 2)), std::sqrt(double(2.)));
  BOOST_CHECK_EQUAL(compute_linear_projection<double>(f, {3, 2, 5, 9}), double(3.));
  F ff = factorize_below(f);
  BOOST_CHECK(ff == F({1, 0}));
  BOOST_CHECK(ff <= f);
  ff = factorize_above(f);
  BOOST_CHECK_EQUAL(ff.num_generators(), 2);
  BOOST_CHECK_EQUAL(ff(0, 0), inf);
  BOOST_CHECK_EQUAL(ff(0, 1), 0);
  BOOST_CHECK_EQUAL(ff(1, 0), 2);
  BOOST_CHECK_EQUAL(ff(1, 1), 1);
  BOOST_CHECK(ff >= f);

  f.add_guaranteed_generator({0, 2});
  BOOST_CHECK_EQUAL(f.num_generators(), 3);
  BOOST_CHECK_EQUAL(f.num_parameters(), 2);

  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(10))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, F({4, 5, 3}, 2)), static_cast<T>(std::sqrt(T(2))));
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {3, 2, 5, 9}), 3);
  ff = factorize_below(f);
  BOOST_CHECK(ff == F({0, 0}));
  BOOST_CHECK(ff <= f);
  ff = factorize_above(f);
  BOOST_CHECK_EQUAL(ff.num_generators(), 3);
  BOOST_CHECK_EQUAL(ff(0, 0), inf);
  BOOST_CHECK_EQUAL(ff(0, 1), 0);
  BOOST_CHECK_EQUAL(ff(1, 0), inf);
  BOOST_CHECK_EQUAL(ff(1, 1), 1);
  BOOST_CHECK_EQUAL(ff(2, 0), 2);
  BOOST_CHECK_EQUAL(ff(2, 1), 2);
  BOOST_CHECK(ff >= f);

  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    T nan = std::numeric_limits<T>::quiet_NaN();
    std::vector<T> v = {nan, 0, 2, 1, nan, 2};
    F f2(v.begin(), v.end(), 2);

    BOOST_CHECK(detail::_is_nan(compute_norm(f2)));
    BOOST_CHECK(detail::_is_nan(compute_euclidean_distance_to(f2, {2, 0})));
    BOOST_CHECK(detail::_is_nan(compute_linear_projection(f2, {3, 0})));
    F f2f = factorize_below(f2);
    BOOST_CHECK_EQUAL(f2f(0, 0), 2);
    BOOST_CHECK_EQUAL(f2f(0, 1), 0);
    f2f = factorize_above(f2);
    BOOST_CHECK_EQUAL(f2f.num_generators(), 3);
    BOOST_CHECK_EQUAL(f2f(0, 0), inf);
    BOOST_CHECK_EQUAL(f2f(0, 1), 0);
    BOOST_CHECK_EQUAL(f2f(1, 0), inf);
    BOOST_CHECK_EQUAL(f2f(1, 1), 1);
    BOOST_CHECK_EQUAL(f2f(2, 0), 2);
    BOOST_CHECK_EQUAL(f2f(2, 1), 2);
  }

  f(0, 0) = 1;
  f(1, 0) = 7;
  f(2, 0) = 5;

  std::vector<std::vector<int>> grid = {{0, 3, 6, 9}, {0, 1, 2, 3}, {0, 4, 8, 16}};
  auto res = compute_coordinates_in_grid(f, grid);
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res.num_generators(), 3);
  BOOST_CHECK_EQUAL(f.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f.num_generators(), 3);
  BOOST_CHECK_EQUAL(res(0, 0), 0);
  BOOST_CHECK_EQUAL(res(1, 0), 2);
  BOOST_CHECK_EQUAL(res(2, 0), 2);

  res = evaluate_coordinates_in_grid(res, grid);
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res.num_generators(), 3);
  BOOST_CHECK_EQUAL(res(0, 0), 0);
  BOOST_CHECK_EQUAL(res(1, 0), 6);
  BOOST_CHECK_EQUAL(res(2, 0), 6);
}

template <class F, typename T>
void test_friends_shifted() {
  T inf = F::has_negative_cones() ? F::T_m_inf : F::T_inf;
  F f({1, 2}, 2);
  f.get_underlying_policy().set_mapping(1, -2);

  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(7))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, F({4, 5, 3}, 2)), static_cast<T>(std::sqrt(T(5))));
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {3, 2, 5, 9}), 4);
  BOOST_CHECK_EQUAL(compute_norm<double>(f), std::sqrt(double(7.)));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to<double>(f, F({4, 5, 3}, 2)), std::sqrt(double(5.)));
  BOOST_CHECK_EQUAL(compute_linear_projection<double>(f, {3, 2, 5, 9}), double(4.));
  F ff = factorize_below(f);
  BOOST_CHECK_EQUAL(ff.num_generators(), 2);
  BOOST_CHECK_EQUAL(ff(0, 0), inf);
  BOOST_CHECK_EQUAL(ff(0, 1), 1);
  BOOST_CHECK_EQUAL(ff(1, 0), 1);
  BOOST_CHECK_EQUAL(ff(1, 1), -1);
  ff = factorize_above(f);
  BOOST_CHECK_EQUAL(ff.num_generators(), 1);
  BOOST_CHECK_EQUAL(ff(0, 0), 2);
  BOOST_CHECK_EQUAL(ff(0, 1), 1);

  f.add_guaranteed_generator({0, -3});
  BOOST_CHECK_EQUAL(f.num_generators(), 3);
  BOOST_CHECK_EQUAL(f.num_parameters(), 2);

  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(16))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, F({4, 5, 3}, 2)), static_cast<T>(std::sqrt(T(5))));
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {3, 2, 5, 9}), -6);
  ff = factorize_below(f);
  BOOST_CHECK_EQUAL(ff.num_generators(), 3);
  BOOST_CHECK_EQUAL(ff(0, 0), inf);
  BOOST_CHECK_EQUAL(ff(0, 1), 1);
  BOOST_CHECK_EQUAL(ff(1, 0), inf);
  BOOST_CHECK_EQUAL(ff(1, 1), -1);
  BOOST_CHECK_EQUAL(ff(2, 0), 0);
  BOOST_CHECK_EQUAL(ff(2, 1), -3);
  ff = factorize_above(f);
  BOOST_CHECK_EQUAL(ff.num_generators(), 1);
  BOOST_CHECK_EQUAL(ff(0, 0), 2);
  BOOST_CHECK_EQUAL(ff(0, 1), 1);

  if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
    T nan = std::numeric_limits<T>::quiet_NaN();
    std::vector<T> v = {nan, 0, 2, 1, nan, 2};
    F f2(v.begin(), v.end(), 2);
    f2.get_underlying_policy().set_mapping(1, -2);

    BOOST_CHECK(detail::_is_nan(compute_norm(f2)));
    BOOST_CHECK(detail::_is_nan(compute_euclidean_distance_to(f2, {2, 0})));
    BOOST_CHECK(detail::_is_nan(compute_linear_projection(f2, {3, 0})));
    F f2f = factorize_below(f2);
    BOOST_CHECK_EQUAL(f2f.num_generators(), 3);
    BOOST_CHECK_EQUAL(f2f(0, 0), inf);
    BOOST_CHECK_EQUAL(f2f(0, 1), 1);
    BOOST_CHECK_EQUAL(f2f(1, 0), inf);
    BOOST_CHECK_EQUAL(f2f(1, 1), -1);
    BOOST_CHECK_EQUAL(f2f(2, 0), 2);
    BOOST_CHECK_EQUAL(f2f(2, 1), -3);
    f2f = factorize_above(f2);
    BOOST_CHECK_EQUAL(f2f.num_generators(), 1);
    BOOST_CHECK_EQUAL(f2f(0, 0), 2);
    BOOST_CHECK_EQUAL(f2f(0, 1), 1);
  }

  f(0, 0) = 1;
  f(1, 0) = 7;
  f(2, 0) = 5;
  f.get_underlying_policy().set_mapping(1, 1);

  std::vector<std::vector<int>> grid = {{0, 3, 6, 9}, {0, 1, 2, 3}, {0, 4, 8, 16}};
  auto res = compute_coordinates_in_grid(f, grid);
  std::cout << res << "\n";
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res.num_generators(), 3);
  BOOST_CHECK_EQUAL(f.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f.num_generators(), 3);
  BOOST_CHECK_EQUAL(res(0, 0), 0);
  BOOST_CHECK_EQUAL(res(1, 0), 2);
  BOOST_CHECK_EQUAL(res(2, 0), 2);
  BOOST_CHECK_EQUAL(res(0, 1), 1);
  BOOST_CHECK_EQUAL(res(1, 1), 2);
  BOOST_CHECK_EQUAL(res(2, 1), 3);

  res = evaluate_coordinates_in_grid(res, grid);
  std::cout << res << "\n";
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res.num_generators(), 3);
  BOOST_CHECK_EQUAL(res(0, 0), 0);
  BOOST_CHECK_EQUAL(res(1, 0), 6);
  BOOST_CHECK_EQUAL(res(2, 0), 6);
  BOOST_CHECK_EQUAL(res(0, 1), 1);
  BOOST_CHECK_EQUAL(res(1, 1), 2);
  BOOST_CHECK_EQUAL(res(2, 1), 3);
}

template <class F, typename T>
void test_friends_1_critical() {
  F f({1, 0});

  BOOST_CHECK_EQUAL(compute_norm(f), static_cast<T>(std::sqrt(T(1))));
  BOOST_CHECK_EQUAL(compute_euclidean_distance_to(f, F({3, 0})), 2);
  BOOST_CHECK_EQUAL(compute_linear_projection(f, {3, 2, 5, 9}), 3);
  BOOST_CHECK(factorize_below(f) == f);
  BOOST_CHECK(factorize_above(f) == f);

  f(0, 0) = 7;

  std::vector<std::vector<int>> grid = {{0, 3, 6, 9}, {0, 1, 2, 3}, {0, 4, 8, 16}};
  auto res = compute_coordinates_in_grid(f, grid);
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(f.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res.num_generators(), 1);
  BOOST_CHECK_EQUAL(f.num_generators(), 1);
  BOOST_CHECK_EQUAL(res(0, 0), 2);
  BOOST_CHECK_EQUAL(res(0, 1), 0);

  res = evaluate_coordinates_in_grid(res, grid);
  BOOST_CHECK_EQUAL(res.num_generators(), 1);
  BOOST_CHECK_EQUAL(res.num_parameters(), 2);
  BOOST_CHECK_EQUAL(res(0, 0), 6);
  BOOST_CHECK_EQUAL(res(0, 1), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_friends, T, list_of_tested_variants) {
  test_friends<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>();
  test_friends_shifted<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>();
  test_friends_1_critical<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>();
}

template <class F, typename T>
void test_unify_intersect(T shift, T step) {
  const int num_param = 2;

  auto get_value = [shift, step](T v) -> T { return shift + (v * step); };

  std::vector<T> v1 = {4, 0, 3, 1, 1, 2, 2, 3};
  F f1(v1.begin(), v1.end(), num_param);
  f1.get_underlying_policy().set_mapping(shift, step);

  std::vector<T> v2 = {5, 0, 2, 1, 1, 2};
  F f2(v2.begin(), v2.end(), num_param);
  f2.get_underlying_policy().set_mapping(shift, step);

  bool modified = unify_lifetimes(f1, f2);
  BOOST_CHECK(modified);
  BOOST_CHECK(f1.num_parameters() == num_param);
  BOOST_CHECK(f1.num_generators() == 4);
  BOOST_CHECK_EQUAL(f1(0, 0), 4);
  BOOST_CHECK_EQUAL(f1(0, 1), get_value(0));
  BOOST_CHECK_EQUAL(f1(1, 0), 2);
  BOOST_CHECK_EQUAL(f1(1, 1), get_value(1));
  BOOST_CHECK_EQUAL(f1(2, 0), 1);
  BOOST_CHECK_EQUAL(f1(2, 1), get_value(2));
  BOOST_CHECK_EQUAL(f1(3, 0), 2);
  BOOST_CHECK_EQUAL(f1(3, 1), get_value(3));

  std::vector<T> v3 = {4, 0, 3, 1, 1, 2, 2, 3};
  F f3(v3.begin(), v3.end(), num_param);
  f3.get_underlying_policy().set_mapping(shift, step);

  modified = intersect_lifetimes(f3, f2);
  BOOST_CHECK(modified);
  BOOST_CHECK(f3.num_parameters() == num_param);
  BOOST_CHECK(f3.num_generators() == 4);
  BOOST_CHECK_EQUAL(f3(0, 0), 5);
  BOOST_CHECK_EQUAL(f3(0, 1), get_value(0));
  BOOST_CHECK_EQUAL(f3(1, 0), 3);
  BOOST_CHECK_EQUAL(f3(1, 1), get_value(1));
  BOOST_CHECK_EQUAL(f3(2, 0), 1);
  BOOST_CHECK_EQUAL(f3(2, 1), get_value(2));
  BOOST_CHECK_EQUAL(f3(3, 0), 1);
  BOOST_CHECK_EQUAL(f3(3, 1), get_value(3));
}

template <class F, typename T>
void test_unify_intersect_1_critical() {
  const int num_param = 2;

  std::vector<T> v1 = {5, 0};
  F f1(v1.begin(), v1.end(), num_param);

  std::vector<T> v2 = {8, 0};
  F f2(v2.begin(), v2.end(), num_param);

  bool modified = unify_lifetimes(f1, f2);
  BOOST_CHECK(!modified);
  BOOST_CHECK(f1.num_parameters() == num_param);
  BOOST_CHECK(f1.num_generators() == 1);
  BOOST_CHECK_EQUAL(f1(0, 0), 5);
  BOOST_CHECK_EQUAL(f1(0, 1), 0);

  modified = intersect_lifetimes(f1, f2);
  BOOST_CHECK(modified);
  BOOST_CHECK(f1.num_parameters() == num_param);
  BOOST_CHECK(f1.num_generators() == 1);
  BOOST_CHECK_EQUAL(f1(0, 0), 8);
  BOOST_CHECK_EQUAL(f1(0, 1), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_unify_intersect, T, list_of_tested_variants) {
  test_unify_intersect<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(0, 1);
  test_unify_intersect<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(-2, 2);
  test_unify_intersect<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(9, -2);
  test_unify_intersect_1_critical<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T>();
}

template <class F, typename T>
void test_serialize(T shift, T step) {
  auto get_value = [shift, step](T v) -> T { return shift + (v * step); };

  std::vector<T> v = {5, 0, 3, 1, 2, 2};
  F f(v.begin(), v.end(), 2);
  f.get_underlying_policy().set_mapping(shift, step);
  BOOST_CHECK(f.num_parameters() == 2);
  BOOST_CHECK(f.num_generators() == 3);
  BOOST_CHECK_EQUAL(f(0, 0), 5);
  BOOST_CHECK_EQUAL(f(0, 1), get_value(0));
  BOOST_CHECK_EQUAL(f(1, 0), 3);
  BOOST_CHECK_EQUAL(f(1, 1), get_value(1));
  BOOST_CHECK_EQUAL(f(2, 0), 2);
  BOOST_CHECK_EQUAL(f(2, 1), get_value(2));

  char* buffer = new char[256];
  std::size_t serializationSize = get_serialization_size_of(f);

  char* ptr = buffer;
  ptr = serialize_value_to_char_buffer(f, ptr);
  BOOST_CHECK_EQUAL(static_cast<std::size_t>(ptr - buffer), serializationSize);

  const char* c_ptr = buffer;
  F f3;
  c_ptr = deserialize_value_from_char_buffer(f3, c_ptr);
  BOOST_CHECK_EQUAL(static_cast<std::size_t>(c_ptr - buffer), serializationSize);
  BOOST_CHECK(f3.num_parameters() == 2);
  BOOST_CHECK(f3.num_generators() == 3);
  BOOST_CHECK_EQUAL(f3(0, 0), 5);
  BOOST_CHECK_EQUAL(f3(0, 1), get_value(0));
  BOOST_CHECK_EQUAL(f3(1, 0), 3);
  BOOST_CHECK_EQUAL(f3(1, 1), get_value(1));
  BOOST_CHECK_EQUAL(f3(2, 0), 2);
  BOOST_CHECK_EQUAL(f3(2, 1), get_value(2));

  delete[] buffer;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_serialize, T, list_of_tested_variants) {
  test_serialize<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(0, 1);
  test_serialize<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(-2, 2);
  test_serialize<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T>(9, -2);
}

template <class F, typename T>
void test_co() {
  F f;
  BOOST_CHECK(f.num_parameters() == 2);
  BOOST_CHECK(f.num_generators() == 1);
  BOOST_CHECK_EQUAL(f(0, 0), F::T_inf);
  BOOST_CHECK_EQUAL(f(0, 1), 0);

  BOOST_CHECK(!f.is_plus_inf());
  BOOST_CHECK(!f.is_minus_inf());
  BOOST_CHECK(!f.is_nan());
  BOOST_CHECK(f.is_finite());

  F f6 = F::minus_inf(2);
  bool change = f6.add_generator({F::T_inf, 0});
  BOOST_CHECK(change);
  BOOST_CHECK_EQUAL(f6(0, 0), F::T_inf);
  BOOST_CHECK_EQUAL(f6(0, 1), 0);

  if constexpr (F::ensures_1_criticality()) {
    std::vector<T> v = {1, 0};
    F f2(v.begin(), v.end(), 2);
    BOOST_CHECK_EQUAL(compute_linear_projection(f2, {3, 2, 5, 9}), 3);
  } else {
    std::vector<T> v = {1, 0, 2, 1, 4, 2};
    F f2(v.begin(), v.end(), 2);
    BOOST_CHECK_EQUAL(compute_linear_projection(f2, {3, 2, 5, 9}), 16);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_co, T, list_of_tested_variants) {
  test_co<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T>();
  test_co<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T>();
}

template <class F, typename T, bool Co>
void test_numerical_limits() {
  const int num_param = 2;

  if constexpr (Co)
    BOOST_CHECK(!std::numeric_limits<F>::has_infinity);
  else
    BOOST_CHECK(std::numeric_limits<F>::has_infinity);
  BOOST_CHECK(std::numeric_limits<F>::has_quiet_NaN);

  BOOST_CHECK(std::numeric_limits<F>::quiet_NaN(num_param).is_nan());
  BOOST_CHECK(std::numeric_limits<F>::minus_infinity(num_param).is_minus_inf());
  if constexpr (Co) {
    BOOST_CHECK_THROW(std::numeric_limits<F>::max(num_param), std::logic_error);
  } else {
    BOOST_CHECK(std::numeric_limits<F>::infinity(num_param).is_plus_inf());
    auto max = std::numeric_limits<F>::max(num_param);
    BOOST_CHECK_EQUAL(max(0, 1), 0);
    BOOST_CHECK_EQUAL(max(0, 0), std::numeric_limits<T>::max());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_numerical_limits, T, list_of_tested_variants) {
  test_numerical_limits<Multi_parameter_filtration_value<Degree_bifiltration<T>>, T, false>();
  test_numerical_limits<Multi_parameter_filtration_value<Degree_bifiltration<T>, false, true>, T, false>();
  test_numerical_limits<Multi_parameter_filtration_value<Degree_bifiltration<T>, true>, T, true>();
  test_numerical_limits<Multi_parameter_filtration_value<Degree_bifiltration<T>, true, true>, T, true>();
}

template <typename T, typename T_alt>
void test_conversions(T shift, T step) {
  auto get_value = [shift, step](T v) -> T { return shift + (v * step); };

  auto test_simplified_value = [&get_value](const auto& f) {
    BOOST_CHECK(f.num_parameters() == 2);
    BOOST_CHECK(f.num_generators() == 2);
    BOOST_CHECK_EQUAL(f(0, 0), 3);
    BOOST_CHECK_EQUAL(f(0, 1), get_value(2));
    BOOST_CHECK_EQUAL(f(1, 0), 5);
    BOOST_CHECK_EQUAL(f(1, 1), get_value(0));
  };

  auto test_value1 = [&get_value](const auto& f) {
    BOOST_CHECK_EQUAL(f.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f.num_generators(), 4);
    BOOST_CHECK_EQUAL(f(0, 0), 5);
    BOOST_CHECK_EQUAL(f(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f(1, 0), 6);
    BOOST_CHECK_EQUAL(f(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f(2, 0), 3);
    BOOST_CHECK_EQUAL(f(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f(3, 0), 4);
    BOOST_CHECK_EQUAL(f(3, 1), get_value(3));
  };

  auto test_value2 = [&get_value](const auto& f) {
    BOOST_CHECK_EQUAL(f.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f.num_generators(), get_value(2) + 1);
    for (std::size_t g = 0; g < static_cast<std::size_t>(get_value(2)) + 1; ++g) {
      BOOST_CHECK_EQUAL(f(g, 1), g);
      if (g == static_cast<std::size_t>(get_value(0)))
        BOOST_CHECK_EQUAL(f(g, 0), 5);
      else if (g == static_cast<std::size_t>(get_value(2)))
        BOOST_CHECK_EQUAL(f(g, 0), 3);
      else
        BOOST_CHECK_EQUAL(f(g, 0), std::decay_t<decltype(f)>::T_inf);
    }
  };

  std::vector<T> v = {5, 6, 3, 4};
  Multi_parameter_filtration_value<Degree_bifiltration<T>> f0(std::move(v), 2);
  f0.get_underlying_policy().set_mapping(shift, step);
  test_value1(f0);
  test_value1(f0.template as_type<T_alt>());

  Multi_parameter_filtration_value<Flat_array_filtration<T>> f1 = f0.template as_type<Flat_array_filtration<T>>();
  test_simplified_value(f1);
  test_value2(f1.template as_type<Degree_bifiltration<T>>());
  test_value2(f1.template as_type<Degree_bifiltration<T_alt>>());

  Multi_parameter_filtration_value<Nested_array_filtration<T>> f2 = f0.template as_type<Nested_array_filtration<T>>();
  test_simplified_value(f2);
  test_value2(f2.template as_type<Degree_bifiltration<T>>());
  test_value2(f2.template as_type<Degree_bifiltration<T_alt>>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_converters, T, list_of_tested_variants) {
  test_conversions<T, T>(0, 1);
  if constexpr (std::is_floating_point_v<T>) {
    test_conversions<T, int>(0, 1);
  } else {
    test_conversions<T, double>(0, 1);
  }

  test_conversions<T, T>(3, 1);
  if constexpr (std::is_floating_point_v<T>) {
    test_conversions<T, int>(3, 1);
  } else {
    test_conversions<T, double>(3, 1);
  }

  test_conversions<T, T>(3, 3);
  if constexpr (std::is_floating_point_v<T>) {
    test_conversions<T, int>(3, 3);
  } else {
    test_conversions<T, double>(3, 3);
  }
}

template <typename T, typename T_alt>
void test_conversions_co(T shift, T step) {
  auto get_value = [shift, step](T v) -> T { return shift + (v * step); };

  auto test_simplified_value = [&get_value](const auto& f) {
    BOOST_CHECK(f.num_parameters() == 2);
    BOOST_CHECK(f.num_generators() == 2);
    BOOST_CHECK_EQUAL(f(0, 0), 4);
    BOOST_CHECK_EQUAL(f(0, 1), get_value(3));
    BOOST_CHECK_EQUAL(f(1, 0), 6);
    BOOST_CHECK_EQUAL(f(1, 1), get_value(1));
  };

  auto test_value1 = [&get_value](const auto& f) {
    BOOST_CHECK(f.num_parameters() == 2);
    BOOST_CHECK(f.num_generators() == 4);
    BOOST_CHECK_EQUAL(f(0, 0), 5);
    BOOST_CHECK_EQUAL(f(0, 1), get_value(0));
    BOOST_CHECK_EQUAL(f(1, 0), 6);
    BOOST_CHECK_EQUAL(f(1, 1), get_value(1));
    BOOST_CHECK_EQUAL(f(2, 0), 3);
    BOOST_CHECK_EQUAL(f(2, 1), get_value(2));
    BOOST_CHECK_EQUAL(f(3, 0), 4);
    BOOST_CHECK_EQUAL(f(3, 1), get_value(3));
  };

  auto test_value2 = [&get_value](const auto& f) {
    BOOST_CHECK_EQUAL(f.num_parameters(), 2);
    BOOST_CHECK_EQUAL(f.num_generators(), get_value(3) + 1);
    for (std::size_t g = 0; g < static_cast<std::size_t>(get_value(3)) + 1; ++g) {
      BOOST_CHECK_EQUAL(f(g, 1), g);
      if (g == static_cast<std::size_t>(get_value(1)))
        BOOST_CHECK_EQUAL(f(g, 0), 6);
      else if (g == static_cast<std::size_t>(get_value(3)))
        BOOST_CHECK_EQUAL(f(g, 0), 4);
      else
        BOOST_CHECK_EQUAL(f(g, 0), std::decay_t<decltype(f)>::T_m_inf);
    }
  };

  std::vector<T> v = {5, 6, 3, 4};
  Multi_parameter_filtration_value<Degree_bifiltration<T>, true> f0(std::move(v), 2);
  f0.get_underlying_policy().set_mapping(shift, step);
  test_value1(f0);
  test_value1(f0.template as_type<T_alt>());

  Multi_parameter_filtration_value<Flat_array_filtration<T>, true> f1 = f0.template as_type<Flat_array_filtration<T>>();
  test_simplified_value(f1);
  test_value2(f1.template as_type<Degree_bifiltration<T>>());
  test_value2(f1.template as_type<Degree_bifiltration<T_alt>>());

  Multi_parameter_filtration_value<Nested_array_filtration<T>, true> f2 =
      f0.template as_type<Nested_array_filtration<T>>();
  test_simplified_value(f2);
  test_value2(f2.template as_type<Degree_bifiltration<T>>());
  test_value2(f2.template as_type<Degree_bifiltration<T_alt>>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(degree_rips_bifiltration_converters_co, T, list_of_tested_variants) {
  test_conversions_co<T, T>(0, 1);
  if constexpr (std::is_floating_point_v<T>) {
    test_conversions_co<T, int>(0, 1);
  } else {
    test_conversions_co<T, double>(0, 1);
  }

  test_conversions_co<T, T>(3, 1);
  if constexpr (std::is_floating_point_v<T>) {
    test_conversions_co<T, int>(3, 1);
  } else {
    test_conversions_co<T, double>(3, 1);
  }

  test_conversions_co<T, T>(3, 3);
  if constexpr (std::is_floating_point_v<T>) {
    test_conversions_co<T, int>(3, 3);
  } else {
    test_conversions_co<T, double>(3, 3);
  }
}
