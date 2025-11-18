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
#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>
#include <gudhi/Degree_rips_bifiltration.h>

using Gudhi::Simplex_tree;
using Gudhi::multi_filtration::Degree_rips_bifiltration;
using Gudhi::multi_filtration::Dynamic_multi_parameter_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration;

template <typename MultiFiltrationValue>
struct Simplex_tree_options_multidimensional_filtration : Gudhi::Simplex_tree_options_default {
  using Filtration_value = MultiFiltrationValue;
};

template <class F>
using Opt = Simplex_tree_options_multidimensional_filtration<F>;

typedef boost::mpl::list<double, float, int> list_of_tested_variants;

// just a very simple test to see if the basics work with the simplex tree.

template <class ST>
void test_multi_simplex_tree()
{
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

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_critical_filtration_io_operator, T, list_of_tested_variants)
{
  std::clog << "Multi_parameter_filtration\n";
  test_multi_simplex_tree<Simplex_tree<Opt<Multi_parameter_filtration<T> > > >();
  std::clog << "Dynamic_multi_parameter_filtration\n";
  test_multi_simplex_tree<Simplex_tree<Opt<Dynamic_multi_parameter_filtration<T> > > >();
  std::clog << "Degree_rips_bifiltration\n";
  test_multi_simplex_tree<Simplex_tree<Opt<Degree_rips_bifiltration<T> > > >();
}
