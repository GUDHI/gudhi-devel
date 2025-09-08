/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <array>
#include <initializer_list>
#include <limits>
#include <utility>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Slicer.h>
#include <gudhi/Thread_safe_slicer.h>
#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>
#include <gudhi/Multi_persistence/Persistence_interface_matrix.h>
#include <gudhi/Multi_persistence/Persistence_interface_cohomology.h>
#include <gudhi/Multi_persistence/Line.h>

using Gudhi::multi_filtration::Dynamic_multi_parameter_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration;
using Gudhi::multi_persistence::Line;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;
using Gudhi::multi_persistence::Persistence_interface_cohomology;
using Gudhi::multi_persistence::Persistence_interface_matrix;
using Gudhi::multi_persistence::Slicer;
using Gudhi::multi_persistence::Thread_safe_slicer;

using T = double;
using Bar = Gudhi::persistence_matrix::Persistence_interval<int, T>;
using Barcode = std::vector<Bar>;
using Flat_barcode = std::vector<std::array<T, 2>>;
using Multi_dimensional_barcode = std::vector<Barcode>;
using Multi_dimensional_flat_barcode = std::vector<Flat_barcode>;
using Test_barcode = std::vector<std::pair<T, T>>;
using Test_multi_dimensional_barcode = std::vector<Test_barcode>;

struct Multi_persistence_r_options
    : Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::INTRUSIVE_SET, true> {
  using Index = std::uint32_t;
  static const bool has_column_pairings = true;
};

struct Multi_persistence_ru_options : Multi_persistence_r_options {
  static const bool has_vine_update = true;
  static const bool can_retrieve_representative_cycles = true;
};

struct Multi_persistence_chain_options : Multi_persistence_ru_options {
  static const bool is_of_boundary_type = false;
  static const Gudhi::persistence_matrix::Column_indexation_types column_indexation_type =
      Gudhi::persistence_matrix::Column_indexation_types::POSITION;
};

typedef boost::mpl::list<
    Slicer<Multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_r_options>>,
    Slicer<Multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_ru_options>>,
    Slicer<Multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_chain_options>>,
    Slicer<Multi_parameter_filtration<T>, Persistence_interface_cohomology<Multi_parameter_filtration<T>>>,
    Slicer<Dynamic_multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_r_options>>,
    Slicer<Dynamic_multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_ru_options>>,
    Slicer<Dynamic_multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_chain_options>>,
    Slicer<Dynamic_multi_parameter_filtration<T>,
           Persistence_interface_cohomology<Dynamic_multi_parameter_filtration<T>>>>
    list_of_tested_variants;

template <class Fil>
Multi_parameter_filtered_complex<Fil> build_rep_cycle_input_complex()
{
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<T>;

  BC bc = {{}, {}, {}, {}, {}, {0, 4}, {1, 4}, {0, 1}, {2, 4}, {3, 4}, {2, 3}, {1, 3}, {5, 6, 7}, {8, 9, 10}};

  DC dc = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2};

  FC fc = {ini{0, 0},
           ini{1, 1},
           ini{2, 2},
           ini{3, 3},
           ini{4, 4},
           ini{5, 5},
           ini{6, 6},
           ini{7, 7},
           ini{8, 8},
           ini{9, 9},
           ini{10, 10},
           ini{11, 11},
           ini{12, 12},
           ini{13, 13}};

  return Complex(bc, dc, fc);
}

template <class Fil>
Multi_parameter_filtered_complex<Fil> build_simple_input_complex()
{
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<T>;

  BC bc = {{}, {}, {}, {0, 1}, {1, 2}, {0, 2}, {3, 4, 5}, {}, {1, 7}};
  DC dc = {0, 0, 0, 1, 1, 1, 2, 0, 1};
  FC fc = {ini{0, 2, 2},
           ini{0, 2, 1},
           ini{0, 1, 3},
           ini{3, 2, 3},
           ini{3, 4, 5},
           ini{6, 3, 5},
           ini{6, 5, 6},
           ini{5, 6, 8},
           ini{5, 7, 8}};

  return Complex(bc, dc, fc);
}

template <class Slicer>
void test_slicer_constructors(const Slicer& s)
{
  BOOST_CHECK_EQUAL(s.get_number_of_cycle_generators(), 9);
  BOOST_CHECK_EQUAL(s.get_number_of_parameters(), 3);
  BOOST_CHECK_EQUAL(s.get_current_order().size(), 9);
  BOOST_CHECK_EQUAL(s.get_slice().size(), 9);
  BOOST_CHECK(!s.get_persistence_algorithm().is_initialized());
}

template <class Slicer>
void test_slicer_constructors_empty(const Slicer& s)
{
  BOOST_CHECK_EQUAL(s.get_number_of_cycle_generators(), 0);
  BOOST_CHECK_EQUAL(s.get_number_of_parameters(), 0);
  BOOST_CHECK(s.get_current_order().empty());
  BOOST_CHECK(s.get_slice().empty());
  BOOST_CHECK(!s.get_persistence_algorithm().is_initialized());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_constructors, Slicer, list_of_tested_variants)
{
  using Fil = typename Slicer::Filtration_value;

  Slicer empty;
  test_slicer_constructors_empty(empty);
  test_slicer_constructors_empty(empty.weak_copy());

  auto cpx = build_simple_input_complex<Fil>();

  Slicer s1(cpx);
  test_slicer_constructors(s1);
  test_slicer_constructors(s1.weak_copy());

  Slicer s2(std::move(cpx));
  BOOST_CHECK_EQUAL(cpx.get_number_of_cycle_generators(), 0);
  test_slicer_constructors(s2);
  test_slicer_constructors(s2.weak_copy());

  Slicer copy(s1);
  test_slicer_constructors(copy);
  test_slicer_constructors(copy.weak_copy());

  Slicer move(std::move(s2));
  test_slicer_constructors_empty(s2);
  test_slicer_constructors_empty(s2.weak_copy());
  test_slicer_constructors(move);
  test_slicer_constructors(move.weak_copy());

  Thread_safe_slicer tss(s1);
  test_slicer_constructors(tss);
  Thread_safe_slicer tssCopy(tss);
  test_slicer_constructors(tssCopy);
}

template <class Slicer>
void test_slicer_accessors(const Slicer& s)
{
  using Fil = typename Slicer::Filtration_value;
  using B = std::vector<typename Slicer::Index>;

  auto box = s.get_bounding_box();
  BOOST_CHECK(box.first == Fil({0, 1, 1}));
  BOOST_CHECK(box.second == Fil({6, 7, 8}));
  BOOST_CHECK_EQUAL(s.get_dimension(0), 0);
  BOOST_CHECK_EQUAL(s.get_dimension(1), 0);
  BOOST_CHECK_EQUAL(s.get_dimension(2), 0);
  BOOST_CHECK_EQUAL(s.get_dimension(3), 1);
  BOOST_CHECK_EQUAL(s.get_dimension(4), 1);
  BOOST_CHECK_EQUAL(s.get_dimension(5), 1);
  BOOST_CHECK_EQUAL(s.get_dimension(6), 2);
  BOOST_CHECK_EQUAL(s.get_dimension(7), 0);
  BOOST_CHECK_EQUAL(s.get_dimension(8), 1);
  BOOST_CHECK(s.get_filtration_value(0) == Fil({0, 2, 2}));
  BOOST_CHECK(s.get_filtration_value(1) == Fil({0, 2, 1}));
  BOOST_CHECK(s.get_filtration_value(2) == Fil({0, 1, 3}));
  BOOST_CHECK(s.get_filtration_value(3) == Fil({3, 2, 3}));
  BOOST_CHECK(s.get_filtration_value(4) == Fil({3, 4, 5}));
  BOOST_CHECK(s.get_filtration_value(5) == Fil({6, 3, 5}));
  BOOST_CHECK(s.get_filtration_value(6) == Fil({6, 5, 6}));
  BOOST_CHECK(s.get_filtration_value(7) == Fil({5, 6, 8}));
  BOOST_CHECK(s.get_filtration_value(8) == Fil({5, 7, 8}));
  BOOST_CHECK(s.get_boundary(0) == B{});
  BOOST_CHECK(s.get_boundary(1) == B{});
  BOOST_CHECK(s.get_boundary(2) == B{});
  BOOST_CHECK(s.get_boundary(3) == (B{0, 1}));
  BOOST_CHECK(s.get_boundary(4) == (B{1, 2}));
  BOOST_CHECK(s.get_boundary(5) == (B{0, 2}));
  BOOST_CHECK(s.get_boundary(6) == (B{3, 4, 5}));
  BOOST_CHECK(s.get_boundary(7) == B{});
  BOOST_CHECK(s.get_boundary(8) == (B{1, 7}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_accessors, Slicer, list_of_tested_variants)
{
  using Fil = typename Slicer::Filtration_value;

  auto cpx = build_simple_input_complex<Fil>();

  Slicer s(cpx);
  test_slicer_accessors(s);
  test_slicer_accessors(s.weak_copy());
  test_slicer_accessors(Thread_safe_slicer(s));
}

template <class Slicer>
void test_slicer_slice_modifiers(Slicer& s)
{
  BOOST_CHECK_EQUAL(s.get_slice().size(), 9);
  s.set_slice({0, 0, 0, 2, 3, 3, 5, 5, 5});
  BOOST_CHECK(s.get_slice() == (std::vector<T>{0, 0, 0, 2, 3, 3, 5, 5, 5}));
  s.push_to(Line<T>({0, 1, 1}));
  BOOST_CHECK(s.get_slice() == (std::vector<T>{1, 1, 2, 3, 4, 6, 6, 7, 7}));
  s.push_to({0, 1, 1});
  BOOST_CHECK(s.get_slice() == (std::vector<T>{1, 1, 2, 3, 4, 6, 6, 7, 7}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_modifiers, Slicer, list_of_tested_variants)
{
  using Fil = typename Slicer::Filtration_value;

  auto cpx = build_simple_input_complex<Fil>();

  Slicer s(cpx);
  test_slicer_slice_modifiers(s);
  Thread_safe_slicer wc = s.weak_copy();
  test_slicer_slice_modifiers(wc);
  Thread_safe_slicer tss(s);
  test_slicer_slice_modifiers(tss);

  std::vector<std::vector<T>> grid = {{0, 2, 4, 8}, {0, 3, 6, 9}, {0, 4, 8, 16}};
  s.coarsen_on_grid(grid, false);

  BOOST_CHECK(s.get_filtration_value(0) == Fil({0, 3, 4}));
  BOOST_CHECK(s.get_filtration_value(1) == Fil({0, 3, 4}));
  BOOST_CHECK(s.get_filtration_value(2) == Fil({0, 3, 4}));
  BOOST_CHECK(s.get_filtration_value(3) == Fil({4, 3, 4}));
  BOOST_CHECK(s.get_filtration_value(4) == Fil({4, 6, 8}));
  BOOST_CHECK(s.get_filtration_value(5) == Fil({8, 3, 8}));
  BOOST_CHECK(s.get_filtration_value(6) == Fil({8, 6, 8}));
  BOOST_CHECK(s.get_filtration_value(7) == Fil({8, 6, 8}));
  BOOST_CHECK(s.get_filtration_value(8) == Fil({8, 9, 8}));

  s.coarsen_on_grid(grid, true);

  BOOST_CHECK(s.get_filtration_value(0) == Fil({0, 1, 1}));
  BOOST_CHECK(s.get_filtration_value(1) == Fil({0, 1, 1}));
  BOOST_CHECK(s.get_filtration_value(2) == Fil({0, 1, 1}));
  BOOST_CHECK(s.get_filtration_value(3) == Fil({2, 1, 1}));
  BOOST_CHECK(s.get_filtration_value(4) == Fil({2, 2, 2}));
  BOOST_CHECK(s.get_filtration_value(5) == Fil({3, 1, 2}));
  BOOST_CHECK(s.get_filtration_value(6) == Fil({3, 2, 2}));
  BOOST_CHECK(s.get_filtration_value(7) == Fil({3, 2, 2}));
  BOOST_CHECK(s.get_filtration_value(8) == Fil({3, 3, 2}));

  Slicer s2(cpx);
  s2.prune_above_dimension(0);
  BOOST_CHECK_EQUAL(s2.get_number_of_cycle_generators(), 4);
  BOOST_CHECK_EQUAL(s2.get_number_of_parameters(), 3);
  BOOST_CHECK_EQUAL(s2.get_dimension(0), 0);
  BOOST_CHECK_EQUAL(s2.get_dimension(1), 0);
  BOOST_CHECK_EQUAL(s2.get_dimension(2), 0);
  BOOST_CHECK_EQUAL(s2.get_dimension(3), 0);
  BOOST_CHECK(s2.get_filtration_value(0) == Fil({0, 2, 1}));
  BOOST_CHECK(s2.get_filtration_value(1) == Fil({0, 2, 2}));
  BOOST_CHECK(s2.get_filtration_value(2) == Fil({0, 1, 3}));
  BOOST_CHECK(s2.get_filtration_value(3) == Fil({5, 6, 8}));
  BOOST_CHECK(s2.get_boundary(0).empty());
  BOOST_CHECK(s2.get_boundary(1).empty());
  BOOST_CHECK(s2.get_boundary(2).empty());
  BOOST_CHECK(s2.get_boundary(3).empty());
}

Test_multi_dimensional_barcode get_barcode(const Barcode& barcode)
{
  Test_multi_dimensional_barcode out(3);  // depends on build_simple_input_complex()...
  for (const auto& bar : barcode) {
    out[bar.dim].emplace_back(bar.birth, bar.death);
  }
  for (auto& dgm : out) std::sort(dgm.begin(), dgm.end());
  return out;
}

Test_multi_dimensional_barcode get_barcode(const Multi_dimensional_barcode& barcode)
{
  Test_multi_dimensional_barcode out(barcode.size());
  for (unsigned int i = 0; i < barcode.size(); ++i) {
    out[i].resize(barcode[i].size());
    for (unsigned int j = 0; j < barcode[i].size(); ++j) {
      out[i][j].first = barcode[i][j].birth;
      out[i][j].second = barcode[i][j].death;
      BOOST_CHECK_EQUAL(i, barcode[i][j].dim);
    }
  }
  for (auto& dgm : out) std::sort(dgm.begin(), dgm.end());
  return out;
}

Test_multi_dimensional_barcode get_barcode(const Multi_dimensional_flat_barcode& barcode)
{
  Test_multi_dimensional_barcode out(barcode.size());
  for (unsigned int i = 0; i < barcode.size(); ++i) {
    out[i].resize(barcode[i].size());
    for (unsigned int j = 0; j < out[i].size(); ++j) {
      out[i][j].first = barcode[i][j][0];
      out[i][j].second = barcode[i][j][1];
    }
  }
  for (auto& dgm : out) std::sort(dgm.begin(), dgm.end());
  return out;
}

Test_barcode get_barcode(const Flat_barcode& barcode)
{
  Test_barcode out(barcode.size());
  for (unsigned int i = 0; i < out.size(); ++i) {
    out[i].first = barcode[i][0];
    out[i].second = barcode[i][1];
  }
  std::sort(out.begin(), out.end());
  return out;
}

template <class Barcode>
void test_barcode(const Barcode& barcode)
{
  auto inf = std::numeric_limits<T>::infinity();
  Test_multi_dimensional_barcode bc = get_barcode(barcode);
  BOOST_CHECK_EQUAL(bc.size(), 3);
  BOOST_CHECK_EQUAL(bc[0].size(), 4);
  BOOST_CHECK(bc[0][0] == (std::pair<T, T>(1, 3)));
  BOOST_CHECK(bc[0][1] == (std::pair<T, T>(1, inf)));
  BOOST_CHECK(bc[0][2] == (std::pair<T, T>(2, 4)));
  BOOST_CHECK(bc[0][3] == (std::pair<T, T>(inf, inf)));
  BOOST_CHECK_EQUAL(bc[1].size(), 1);
  BOOST_CHECK(bc[1][0] == (std::pair<T, T>(6, inf)));
  BOOST_CHECK_EQUAL(bc[2].size(), 0);
}

template <>
void test_barcode(const Flat_barcode& barcode)
{
  auto inf = std::numeric_limits<T>::infinity();
  Test_barcode bc = get_barcode(barcode);
  BOOST_CHECK_EQUAL(bc.size(), 5);
  BOOST_CHECK(bc[0] == (std::pair<T, T>(1, 3)));
  BOOST_CHECK(bc[1] == (std::pair<T, T>(1, inf)));
  BOOST_CHECK(bc[2] == (std::pair<T, T>(2, 4)));
  BOOST_CHECK(bc[3] == (std::pair<T, T>(6, inf)));
  BOOST_CHECK(bc[4] == (std::pair<T, T>(inf, inf)));
}

template <class Barcode>
void test_barcode_ignore_inf(const Barcode& barcode)
{
  auto inf = std::numeric_limits<T>::infinity();
  Test_multi_dimensional_barcode bc = get_barcode(barcode);
  BOOST_CHECK_EQUAL(bc.size(), 3);
  BOOST_CHECK_EQUAL(bc[0].size(), 3);
  BOOST_CHECK(bc[0][0] == (std::pair<T, T>(1, 3)));
  BOOST_CHECK(bc[0][1] == (std::pair<T, T>(1, inf)));
  BOOST_CHECK(bc[0][2] == (std::pair<T, T>(2, 4)));
  BOOST_CHECK_EQUAL(bc[1].size(), 1);
  BOOST_CHECK(bc[1][0] == (std::pair<T, T>(6, inf)));
  BOOST_CHECK_EQUAL(bc[2].size(), 0);
}

template <>
void test_barcode_ignore_inf(const Flat_barcode& barcode)
{
  auto inf = std::numeric_limits<T>::infinity();
  Test_barcode bc = get_barcode(barcode);
  BOOST_CHECK_EQUAL(bc.size(), 4);
  BOOST_CHECK(bc[0] == (std::pair<T, T>(1, 3)));
  BOOST_CHECK(bc[1] == (std::pair<T, T>(1, inf)));
  BOOST_CHECK(bc[2] == (std::pair<T, T>(2, 4)));
  BOOST_CHECK(bc[3] == (std::pair<T, T>(6, inf)));
}

template <class Slicer>
void test_slicer_persistence(Slicer& s)
{
  using Index = typename Slicer::Index;
  T inf = std::numeric_limits<T>::infinity();

  s.set_slice({1, 2, 1, 6, 4, 3, inf, inf, inf});
  s.initialize_persistence_computation(false);

  BOOST_CHECK(s.persistence_computation_is_initialized());
  BOOST_CHECK(s.get_current_order() == (std::vector<Index>{0, 2, 1, 7, 5, 4, 3, 8, 6}));

  test_barcode(s.template get_barcode<true>());
  test_barcode(s.template get_barcode<false>());
  test_barcode(s.template get_flat_barcode<true>());
  test_barcode(s.template get_flat_barcode<false>());

  s.initialize_persistence_computation(true);

  BOOST_CHECK(s.persistence_computation_is_initialized());
  BOOST_CHECK(s.get_current_order() == (std::vector<Index>{0, 2, 1, 5, 4, 3}));

  test_barcode_ignore_inf(s.template get_barcode<true>());
  test_barcode_ignore_inf(s.template get_barcode<false>());
  test_barcode_ignore_inf(s.template get_flat_barcode<true>());
  test_barcode_ignore_inf(s.template get_flat_barcode<false>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_persistence, Slicer, list_of_tested_variants)
{
  using Fil = typename Slicer::Filtration_value;

  auto cpx = build_simple_input_complex<Fil>();

  Slicer s(cpx);
  test_slicer_persistence(s);
  Thread_safe_slicer wc = s.weak_copy();
  test_slicer_persistence(wc);
  Thread_safe_slicer tss(s);
  test_slicer_persistence(tss);
}

typedef boost::mpl::list<
    Slicer<Multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_ru_options>>,
    Slicer<Multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_chain_options>>,
    Slicer<Dynamic_multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_ru_options>>,
    Slicer<Dynamic_multi_parameter_filtration<T>, Persistence_interface_matrix<Multi_persistence_chain_options>>>
    list_of_tested_variants2;

template <class Slicer>
void test_slicer_vineyard(Slicer& s)
{
  auto inf = std::numeric_limits<T>::infinity();

  s.set_slice({1, 2, 1, 6, 4, 3, inf, inf, inf});
  s.initialize_persistence_computation(false);

  Test_multi_dimensional_barcode bc = get_barcode(s.get_barcode());
  BOOST_CHECK_EQUAL(bc.size(), 3);
  BOOST_CHECK_EQUAL(bc[0].size(), 4);
  BOOST_CHECK(bc[0][0] == (std::pair<T, T>(1, 3)));
  BOOST_CHECK(bc[0][1] == (std::pair<T, T>(1, inf)));
  BOOST_CHECK(bc[0][2] == (std::pair<T, T>(2, 4)));
  BOOST_CHECK(bc[0][3] == (std::pair<T, T>(inf, inf)));
  BOOST_CHECK_EQUAL(bc[1].size(), 1);
  BOOST_CHECK(bc[1][0] == (std::pair<T, T>(6, inf)));
  BOOST_CHECK_EQUAL(bc[2].size(), 0);

  s.push_to(Line<T>({0, 1, 1}));
  s.vineyard_update();

  bc = get_barcode(s.get_barcode());
  BOOST_CHECK_EQUAL(bc.size(), 3);
  BOOST_CHECK_EQUAL(bc[0].size(), 4);
  BOOST_CHECK(bc[0][0] == (std::pair<T, T>(1, 3)));
  BOOST_CHECK(bc[0][1] == (std::pair<T, T>(1, inf)));
  BOOST_CHECK(bc[0][2] == (std::pair<T, T>(2, 4)));
  BOOST_CHECK(bc[0][3] == (std::pair<T, T>(7, 7)));
  BOOST_CHECK_EQUAL(bc[1].size(), 1);
  BOOST_CHECK(bc[1][0] == (std::pair<T, T>(6, 6)));
  BOOST_CHECK_EQUAL(bc[2].size(), 0);

  s.push_to(Line<T>({2, 4, 2}));
  s.vineyard_update();

  bc = get_barcode(s.get_barcode());
  BOOST_CHECK_EQUAL(bc.size(), 3);
  BOOST_CHECK_EQUAL(bc[0].size(), 4);
  BOOST_CHECK(bc[0][0] == (std::pair<T, T>(-1, inf)));
  BOOST_CHECK(bc[0][1] == (std::pair<T, T>(0, 1)));
  BOOST_CHECK(bc[0][2] == (std::pair<T, T>(1, 3)));
  BOOST_CHECK(bc[0][3] == (std::pair<T, T>(6, 6)));
  BOOST_CHECK_EQUAL(bc[1].size(), 1);
  BOOST_CHECK(bc[1][0] == (std::pair<T, T>(4, 4)));
  BOOST_CHECK_EQUAL(bc[2].size(), 0);

  s.push_to(Line<T>({2, 1, 0}));
  s.vineyard_update();

  bc = get_barcode(s.get_barcode());
  BOOST_CHECK_EQUAL(bc.size(), 3);
  BOOST_CHECK_EQUAL(bc[0].size(), 4);
  BOOST_CHECK(bc[0][0] == (std::pair<T, T>(1, inf)));
  BOOST_CHECK(bc[0][1] == (std::pair<T, T>(2, 3)));
  BOOST_CHECK(bc[0][2] == (std::pair<T, T>(3, 5)));
  BOOST_CHECK(bc[0][3] == (std::pair<T, T>(8, 8)));
  BOOST_CHECK_EQUAL(bc[1].size(), 1);
  BOOST_CHECK(bc[1][0] == (std::pair<T, T>(5, 6)));
  BOOST_CHECK_EQUAL(bc[2].size(), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_vineyard, Slicer, list_of_tested_variants2)
{
  using Fil = typename Slicer::Filtration_value;

  auto cpx = build_simple_input_complex<Fil>();

  Slicer s(cpx);
  test_slicer_vineyard(s);
  Thread_safe_slicer wc = s.weak_copy();
  test_slicer_vineyard(wc);
  Thread_safe_slicer tss(s);
  test_slicer_vineyard(tss);
}

template <class Slicer>
void test_slicer_rep_cycles(Slicer& s)
{
  using Cycle = typename Slicer::Cycle;

  s.set_slice({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13});
  s.initialize_persistence_computation(false);

  std::vector<std::vector<Cycle>> cycles = s.get_representative_cycles();

  BOOST_CHECK_EQUAL(cycles.size(), 3);

  BOOST_CHECK_EQUAL(cycles[0].size(), 5);
  BOOST_CHECK(cycles[0][0] == Cycle({0}));
  BOOST_CHECK(cycles[0][1] == Cycle({0, 1}));
  BOOST_CHECK(cycles[0][2] == Cycle({0, 2}));
  BOOST_CHECK(cycles[0][3] == Cycle({0, 3}));
  BOOST_CHECK(cycles[0][4] == Cycle({0, 4}));

  BOOST_CHECK_EQUAL(cycles[1].size(), 3);
  BOOST_CHECK(cycles[1][0] == Cycle({5, 6, 7}));
  BOOST_CHECK(cycles[1][1] == Cycle({8, 9, 10}));
  BOOST_CHECK(cycles[1][2] == Cycle({6, 9, 11}));

  BOOST_CHECK_EQUAL(cycles[2].size(), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_rep_cycles, Slicer, list_of_tested_variants2)
{
  using Fil = typename Slicer::Filtration_value;

  auto cpx = build_rep_cycle_input_complex<Fil>();

  Slicer s(cpx);
  test_slicer_rep_cycles(s);
  Thread_safe_slicer wc = s.weak_copy();
  test_slicer_rep_cycles(wc);
  Thread_safe_slicer tss(s);
  test_slicer_rep_cycles(tss);
}