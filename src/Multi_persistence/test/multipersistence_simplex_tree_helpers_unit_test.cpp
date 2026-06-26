/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <initializer_list>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/multi_simplex_tree_helpers.h>
#include <gudhi/Multi_parameter_filtration.h>

using Gudhi::multi_filtration::Multi_parameter_filtration;
using Gudhi::multi_persistence::make_multi_dimensional;
using Gudhi::multi_persistence::make_one_dimensional;
using Gudhi::multi_persistence::Simplex_tree_options_multidimensional_filtration;

template <typename T>
using Multi_options = Simplex_tree_options_multidimensional_filtration<Multi_parameter_filtration<T> >;
template <typename T>
using Std_options = Simplex_tree_options_multidimensional_filtration<T>;
template <typename T>
using Multi_tree = Gudhi::Simplex_tree<Multi_options<T> >;
template <typename T>
using Std_tree = Gudhi::Simplex_tree<Std_options<T> >;

template <typename T>
Multi_tree<T> build_multi_simplex_tree()
{
  Multi_tree<T> st;
  st.insert_simplex_and_subfaces({0, 1, 2}, std::initializer_list<T>{0, 1, 2});
  st.insert_simplex_and_subfaces({0, 1, 3}, std::initializer_list<T>{3, 4, 5});
  st.insert_simplex_and_subfaces({2, 3, 4}, std::initializer_list<T>{6, 7, 8});
  st.set_num_parameters(3);
  return st;
}

template <typename T>
Std_tree<T> build_std_simplex_tree()
{
  Std_tree<T> st;
  st.insert_simplex_and_subfaces({0, 1, 2}, 1);
  st.insert_simplex_and_subfaces({0, 1, 3}, 4);
  st.insert_simplex_and_subfaces({2, 3, 4}, 7);
  return st;
}

template <typename T>
Multi_tree<T> build_proj_multi_simplex_tree()
{
  Multi_tree<T> st;
  st.insert_simplex_and_subfaces({0, 1, 2}, std::initializer_list<T>{0, 1, 2});
  st.insert_simplex_and_subfaces({0, 1, 3}, std::initializer_list<T>{0, 4, 2});
  st.insert_simplex_and_subfaces({2, 3, 4}, std::initializer_list<T>{0, 7, 2});
  st.set_num_parameters(3);
  return st;
}

BOOST_AUTO_TEST_CASE(multi_simplex_tree_helpers)
{
  Multi_tree<double> multi_st = build_multi_simplex_tree<double>();
  Std_tree<int> std_st = build_std_simplex_tree<int>();
  Multi_tree<float> proj_multi_st = build_proj_multi_simplex_tree<float>();

  Std_tree<int> test1 = make_one_dimensional<Std_options<int> >(multi_st, 1);
  BOOST_CHECK(test1 == std_st);

  Multi_tree<float> test2 =
      make_multi_dimensional<Multi_options<float> >(std_st, std::initializer_list<float>{0, 3, 2}, 1);
  BOOST_CHECK(test2 == proj_multi_st);
}
