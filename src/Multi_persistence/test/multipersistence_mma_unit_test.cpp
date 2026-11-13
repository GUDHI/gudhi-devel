/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstdint>
#include <initializer_list>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/multiparameter_module_approximation.h>
#include <gudhi/Multi_filtration/Flat_array_filtration.h>
#include <gudhi/Multi_filtration/Nested_array_filtration.h>
#include <gudhi/Multi_parameter_filtration_value.h>
#include <gudhi/Slicer.h>
#include <gudhi/Multi_persistence/Persistence_interface_vineyard.h>

using Gudhi::multi_filtration::Flat_array_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration_value;
using Gudhi::multi_filtration::Nested_array_filtration;
using Gudhi::multi_persistence::Box;
using Gudhi::multi_persistence::Module;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;
using Gudhi::multi_persistence::multiparameter_module_approximation;
using Gudhi::multi_persistence::Persistence_interface_vineyard;
using Gudhi::multi_persistence::Slicer;

struct Multi_persistence_vineyard_ru_options : Gudhi::vineyard::Default_vineyard_options {
  static constexpr bool is_RU = true;
};

using I = std::uint32_t;
using D = int;

using list_of_tested_variants =
    boost::mpl::list<Slicer<Multi_parameter_filtration_value<Flat_array_filtration<double>>,
                            Persistence_interface_vineyard<Multi_persistence_vineyard_ru_options>>,
                     Slicer<Multi_parameter_filtration_value<Nested_array_filtration<double>>,
                            Persistence_interface_vineyard<Multi_persistence_vineyard_ru_options>>,
                     Slicer<Multi_parameter_filtration_value<Flat_array_filtration<float>>,
                            Persistence_interface_vineyard<Multi_persistence_vineyard_ru_options>>,
                     Slicer<Multi_parameter_filtration_value<Nested_array_filtration<float>>,
                            Persistence_interface_vineyard<Multi_persistence_vineyard_ru_options>>>;

template <class Fil>
Multi_parameter_filtered_complex<Fil, I, D> build_simple_input_complex() {
  using Complex = Multi_parameter_filtered_complex<Fil, I, D>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<typename Fil::value_type>;

  BC bc = {{}, {}, {}, {0, 1}, {1, 2}, {0, 2}, {3, 4, 5}, {}, {1, 7}};
  DC dc = {0, 0, 0, 1, 1, 1, 2, 0, 1};
  FC fc = {ini{0, 2, 2}, ini{0, 2, 1}, ini{0, 1, 3}, ini{3, 2, 3}, ini{3, 4, 5},
           ini{6, 3, 5}, ini{6, 5, 6}, ini{5, 6, 8}, ini{5, 7, 8}};

  return Complex(bc, dc, fc);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(mma, Slicer_t, list_of_tested_variants) {
  using Fil = typename Slicer_t::Filtration_value;
  using T = typename Fil::value_type;
  using S = typename Module<T>::Summand_t;

  T inf = S::T_inf;

  const int numParam = 3;
  S s1({0, 2, 1, 1, 1, 3}, {inf, inf, inf}, numParam, 0);
  S s2({0, 2, 2}, {3, inf, inf, inf, inf, 3}, numParam, 0);
  S s3({0, 2, 3}, {3, inf, inf, inf, 4, inf, inf, inf, 5}, numParam, 0);
  S s4({6, 6, 8}, {inf, 7, inf}, numParam, 0);
  S s5({6, 4, 5}, {inf, 5, inf, inf, inf, 6}, numParam, 1);

  auto cpx = build_simple_input_complex<Fil>();
  Slicer_t slicer(cpx);

  auto mod = multiparameter_module_approximation(slicer, 1, {0, 1, 0}, {5, 8, 8});
  BOOST_CHECK_EQUAL(mod.get_number_of_parameters(), numParam);
  BOOST_CHECK_EQUAL(mod.size(), 5);
  BOOST_CHECK(mod.get_summand(0) == s1);
  BOOST_CHECK(mod.get_summand(1) == s2);
  BOOST_CHECK(mod.get_summand(2) == s3);
  BOOST_CHECK(mod.get_summand(3) == s4);
  BOOST_CHECK(mod.get_summand(4) == s5);

  S s6({0, 2, 1, 1, 1, 3}, {5, 8, 8}, numParam, 0);
  S s7({0, 2, 2}, {3, 8, 8, 5, 8, 3}, numParam, 0);
  S s8({0, 2, 3}, {3, 8, 8, 5, 4, 8, 5, 8, 5}, numParam, 0);
  S s9({5, 6, 8}, {5, 7, 8}, numParam, 0);
  S s10({6, 4, 8, 6, 4.5, 6, 6.5, 4.25, 5}, {5, 5, 8, 5, 8, 6}, numParam, 1);

  mod = multiparameter_module_approximation(slicer, 1, {0, 1, 0}, {5, 8, 8}, {2, 1, 4}, true);
  BOOST_CHECK_EQUAL(mod.get_number_of_parameters(), numParam);
  BOOST_CHECK_EQUAL(mod.size(), 5);
  BOOST_CHECK(mod.get_summand(0) == s6);
  BOOST_CHECK(mod.get_summand(1) == s8);
  BOOST_CHECK(mod.get_summand(2) == s7);
  BOOST_CHECK(mod.get_summand(3) == s9);
  BOOST_CHECK(mod.get_summand(4) == s10);

  S s11({0, 2, 1}, {5, 8, 8}, numParam, 0);
  S s12({0, 2, 2}, {5, 8, 8}, numParam, 0);
  S s13({0, 2, 3}, {5, 8, 8}, numParam, 0);
  S s14({5, 6, 8}, {5, 8, 8}, numParam, 0);
  S s15({6, 4, 5}, {5, 8, 8}, numParam, 1);

  mod = multiparameter_module_approximation(slicer, 1, {0, 1, 0}, {5, 8, 8}, {1, 0, 0}, true);
  BOOST_CHECK_EQUAL(mod.get_number_of_parameters(), numParam);
  BOOST_CHECK_EQUAL(mod.size(), 5);
  BOOST_CHECK(mod.get_summand(0) == s11);
  BOOST_CHECK(mod.get_summand(1) == s12);
  BOOST_CHECK(mod.get_summand(2) == s13);
  BOOST_CHECK(mod.get_summand(3) == s14);
  BOOST_CHECK(mod.get_summand(4) == s15);
}
