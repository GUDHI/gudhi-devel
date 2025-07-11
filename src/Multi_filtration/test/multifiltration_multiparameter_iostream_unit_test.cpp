/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <sstream>
#include <type_traits>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_filtration"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>
#include <gudhi/Degree_rips_bifiltration.h>

using Gudhi::multi_filtration::Degree_rips_bifiltration;
using Gudhi::multi_filtration::Dynamic_multi_parameter_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration;

typedef boost::mpl::list<double, float, int> list_of_tested_variants;

template <class F>
void test_io_operator(const F& f)
{
  std::stringstream ss;
  ss << f;
  F f2;
  ss >> f2;
  BOOST_CHECK(f == f2 || (f.is_nan() && f2.is_nan()));
}

template <class F, typename T>
void test_io()
{
  const int num_param = 2;
  std::vector<T> v1, v2 = {0, 0, 1, -2, 2, -1, 3, 0};

  F f1(v1.begin(), v1.end(), num_param);
  F f2(v2.begin(), v2.end(), num_param);
  F f3 = F::inf(num_param);
  F f4 = F::minus_inf(num_param);

  test_io_operator(f1);
  test_io_operator(f2);
  test_io_operator(f3);
  test_io_operator(f4);

  if constexpr (!std::is_same_v<T, int>) {
    F f5 = F::nan(num_param);
    test_io_operator(f5);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_io_operator, T, list_of_tested_variants)
{
  std::clog << "Multi_parameter_filtration\n";
  test_io<Multi_parameter_filtration<T>, T>();
  std::clog << "Dynamic_multi_parameter_filtration\n";
  test_io<Dynamic_multi_parameter_filtration<T>, T>();
  std::clog << "Degree_rips_bifiltration\n";
  test_io<Degree_rips_bifiltration<T>, T>();
}
