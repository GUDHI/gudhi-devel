/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_filtration"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Multi_filtration/Degree_bifiltration.h>
#include <gudhi/Multi_filtration/Flat_array_filtration.h>
#include <gudhi/Multi_filtration/Nested_array_filtration.h>
#include <gudhi/Multi_parameter_filtration_value.h>

using Gudhi::Simplex_tree;
using Gudhi::multi_filtration::Degree_bifiltration;
using Gudhi::multi_filtration::Flat_array_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration_value;
using Gudhi::multi_filtration::Nested_array_filtration;

template <typename MultiFiltrationValue>
struct Simplex_tree_options_multidimensional_filtration : Gudhi::Simplex_tree_options_default {
  using Filtration_value = MultiFiltrationValue;
};

template <class F>
using Opt = Simplex_tree_options_multidimensional_filtration<F>;

using list_of_tested_variants = boost::mpl::list<double, float, int>;

// just a very simple test to see if the basics work with the simplex tree.

template <class ST>
void test_multi_simplex_tree() {
  using F = typename ST::Filtration_value;
  using T = typename F::value_type;
  using ini = std::initializer_list<T>;

  ST simplexTree;

  simplexTree.insert_simplex_and_subfaces({0, 1, 2}, ini{3, 0});
  simplexTree.insert_simplex_and_subfaces({1, 3}, ini{4, 0});
  simplexTree.insert_simplex_and_subfaces({4, 5}, ini{6, 0});
  simplexTree.insert_simplex_and_subfaces({3, 4, 5, 6}, ini{5, 0});
  simplexTree.insert_simplex_and_subfaces({2, 6}, ini{7, 0});
  simplexTree.insert_simplex_and_subfaces({3, 4}, ini{8, 0});
  simplexTree.insert_simplex_and_subfaces({0, 1, 2}, ini{2, 0});
  simplexTree.insert_simplex_and_subfaces({4, 5, 6}, ini{4, 0});
  simplexTree.insert_simplex_and_subfaces({2, 6}, ini{1, 0});
  simplexTree.insert_simplex_and_subfaces({1, 3}, ini{8, 0});

  BOOST_CHECK((simplexTree.filtration(simplexTree.find({0})) == F({2, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({1})) == F({2, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({2})) == F({1, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({3})) == F({4, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({4})) == F({4, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({5})) == F({4, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({6})) == F({1, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({0, 1})) == F({2, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({0, 2})) == F({2, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({1, 2})) == F({2, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({3, 4})) == F({5, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({3, 5})) == F({5, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({3, 6})) == F({5, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({4, 5})) == F({4, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({4, 6})) == F({4, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({5, 6})) == F({4, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({1, 3})) == F({4, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({2, 6})) == F({1, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({0, 1, 2})) == F({2, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({3, 4, 5})) == F({5, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({3, 4, 6})) == F({5, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({3, 5, 6})) == F({5, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({4, 5, 6})) == F({4, 0})));
  BOOST_CHECK((simplexTree.filtration(simplexTree.find({3, 4, 5, 6})) == F({5, 0})));

  std::clog << simplexTree << "\n";
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_io_operator, T, list_of_tested_variants) {
  std::clog << "Flat_array_filtration\n";
  test_multi_simplex_tree<Simplex_tree<Opt<Multi_parameter_filtration_value<Flat_array_filtration<T> > > > >();
  std::clog << "Nested_array_filtration\n";
  test_multi_simplex_tree<Simplex_tree<Opt<Multi_parameter_filtration_value<Nested_array_filtration<T> > > > >();
  std::clog << "Degree_bifiltration\n";
  test_multi_simplex_tree<Simplex_tree<Opt<Multi_parameter_filtration_value<Degree_bifiltration<T> > > > >();
}
