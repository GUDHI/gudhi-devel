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
#include <cstddef>
#include <initializer_list>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/simple_mdspan.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_filtration/Flat_array_filtration.h>
#include <gudhi/Multi_filtration/Nested_array_filtration.h>
#include <gudhi/Multi_parameter_filtration_value.h>
#include <gudhi/slicer_helpers.h>
#include <gudhi/Slicer.h>
#include <gudhi/Thread_safe_slicer.h>
#include <gudhi/multi_simplex_tree_helpers.h>
#include <gudhi/multi_persistence_landscapes.h>
#include <gudhi/Multi_persistence/Persistence_interface_homology.h>
#include <gudhi/Multi_persistence/Persistence_interface_cohomology.h>
#include <gudhi/Multi_persistence/Persistence_interface_vineyard.h>

using Gudhi::multi_filtration::Flat_array_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration_value;
using Gudhi::multi_filtration::Nested_array_filtration;
using Gudhi::multi_persistence::compute_slicer_landscapes_on_grid;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;
using Gudhi::multi_persistence::Persistence_interface_cohomology;
using Gudhi::multi_persistence::Persistence_interface_homology;
using Gudhi::multi_persistence::Persistence_interface_vineyard;
using Gudhi::multi_persistence::persistence_on_slices;
using Gudhi::multi_persistence::Slicer;
using Gudhi::multi_persistence::Thread_safe_slicer;

using I = std::uint32_t;
using D = int;
using T = double;
using Flat_barcode = std::vector<std::array<T, 2>>;
using Multi_dimensional_flat_barcode = std::vector<Flat_barcode>;
using Test_barcode = std::vector<std::pair<T, T>>;
using Test_multi_dimensional_barcode = std::vector<Test_barcode>;

struct Multi_persistence_r_options
    : Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::INTRUSIVE_SET, true> {
  using Index = I;
  using Dimension = D;
  static const bool has_column_pairings = true;
};

struct Multi_persistence_ru_options : Multi_persistence_r_options {
  static const bool can_retrieve_representative_cycles = true;
};

struct Multi_persistence_chain_options : Multi_persistence_ru_options {
  static const bool is_of_boundary_type = false;
  static const Gudhi::persistence_matrix::Column_indexation_types column_indexation_type =
      Gudhi::persistence_matrix::Column_indexation_types::POSITION;
};

struct Multi_persistence_vineyard_ru_options : Gudhi::vineyard::Default_vineyard_options {
  static constexpr bool is_RU = true;
};

struct Multi_persistence_vineyard_chain_options : Gudhi::vineyard::Default_vineyard_options {
  static constexpr bool is_RU = false;
};

using Flat_MFV = Multi_parameter_filtration_value<Flat_array_filtration<T>>;
using Nested_MFV = Multi_parameter_filtration_value<Nested_array_filtration<T>>;

using list_of_tested_slicer_variants =
    boost::mpl::list<Slicer<Flat_MFV, Persistence_interface_homology<Multi_persistence_r_options, Flat_MFV>>,
                     Slicer<Flat_MFV, Persistence_interface_homology<Multi_persistence_ru_options, Flat_MFV>>,
                     Slicer<Flat_MFV, Persistence_interface_homology<Multi_persistence_chain_options, Flat_MFV>>,
                     Slicer<Flat_MFV, Persistence_interface_cohomology<Flat_MFV>>,
                     Slicer<Flat_MFV, Persistence_interface_vineyard<Multi_persistence_vineyard_ru_options>>,
                     Slicer<Flat_MFV, Persistence_interface_vineyard<Multi_persistence_vineyard_chain_options>>,
                     Slicer<Nested_MFV, Persistence_interface_homology<Multi_persistence_r_options, Nested_MFV>>,
                     Slicer<Nested_MFV, Persistence_interface_homology<Multi_persistence_ru_options, Nested_MFV>>,
                     Slicer<Nested_MFV, Persistence_interface_homology<Multi_persistence_chain_options, Nested_MFV>>,
                     Slicer<Nested_MFV, Persistence_interface_cohomology<Nested_MFV>>,
                     Slicer<Nested_MFV, Persistence_interface_vineyard<Multi_persistence_vineyard_ru_options>>,
                     Slicer<Nested_MFV, Persistence_interface_vineyard<Multi_persistence_vineyard_chain_options>>>;

template <class Fil>
Multi_parameter_filtered_complex<Fil, I, D> build_simple_input_complex() {
  using Complex = Multi_parameter_filtered_complex<Fil, I, D>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<T>;

  BC bc = {{}, {}, {}, {0, 1}, {1, 2}, {0, 2}, {3, 4, 5}, {}, {1, 7}};
  DC dc = {0, 0, 0, 1, 1, 1, 2, 0, 1};
  FC fc = {ini{0, 2, 2}, ini{0, 2, 1}, ini{0, 1, 3}, ini{3, 2, 3}, ini{3, 4, 5},
           ini{6, 3, 5}, ini{6, 5, 6}, ini{5, 6, 8}, ini{5, 7, 8}};

  return Complex(bc, dc, fc);
}

Test_multi_dimensional_barcode get_barcode(const Multi_dimensional_flat_barcode& barcode) {
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

template <class Slicer>
void test_slicer_batch_persistence(Slicer& s) {
  using Index = typename Slicer::Index;
  T inf = std::numeric_limits<T>::infinity();

  auto barcodes = persistence_on_slices<T>(s, {{0, 1, 1}, {2, 4, 2}, {2, 1, 0}}, {}, false);
  std::vector<Test_multi_dimensional_barcode> orderedBarcodes(barcodes.size());
  Index i = 0;
  for (const auto& dgm : barcodes) {
    orderedBarcodes[i] = get_barcode(dgm);
    ++i;
  }

  BOOST_CHECK_EQUAL(orderedBarcodes.size(), 3);

  BOOST_CHECK_EQUAL(orderedBarcodes[0].size(), 3);
  BOOST_CHECK_EQUAL(orderedBarcodes[0][0].size(), 4);
  BOOST_CHECK(orderedBarcodes[0][0][0] == (std::pair<T, T>(1, 3)));
  BOOST_CHECK(orderedBarcodes[0][0][1] == (std::pair<T, T>(1, inf)));
  BOOST_CHECK(orderedBarcodes[0][0][2] == (std::pair<T, T>(2, 4)));
  BOOST_CHECK(orderedBarcodes[0][0][3] == (std::pair<T, T>(7, 7)));
  BOOST_CHECK_EQUAL(orderedBarcodes[0][1].size(), 1);
  BOOST_CHECK(orderedBarcodes[0][1][0] == (std::pair<T, T>(6, 6)));
  BOOST_CHECK_EQUAL(orderedBarcodes[0][2].size(), 0);

  BOOST_CHECK_EQUAL(orderedBarcodes[1].size(), 3);
  BOOST_CHECK_EQUAL(orderedBarcodes[1][0].size(), 4);
  BOOST_CHECK(orderedBarcodes[1][0][0] == (std::pair<T, T>(-1, inf)));
  BOOST_CHECK(orderedBarcodes[1][0][1] == (std::pair<T, T>(0, 1)));
  BOOST_CHECK(orderedBarcodes[1][0][2] == (std::pair<T, T>(1, 3)));
  BOOST_CHECK(orderedBarcodes[1][0][3] == (std::pair<T, T>(6, 6)));
  BOOST_CHECK_EQUAL(orderedBarcodes[1][1].size(), 1);
  BOOST_CHECK(orderedBarcodes[1][1][0] == (std::pair<T, T>(4, 4)));
  BOOST_CHECK_EQUAL(orderedBarcodes[1][2].size(), 0);

  BOOST_CHECK_EQUAL(orderedBarcodes[2].size(), 3);
  BOOST_CHECK_EQUAL(orderedBarcodes[2][0].size(), 4);
  BOOST_CHECK(orderedBarcodes[2][0][0] == (std::pair<T, T>(1, inf)));
  BOOST_CHECK(orderedBarcodes[2][0][1] == (std::pair<T, T>(2, 3)));
  BOOST_CHECK(orderedBarcodes[2][0][2] == (std::pair<T, T>(3, 5)));
  BOOST_CHECK(orderedBarcodes[2][0][3] == (std::pair<T, T>(8, 8)));
  BOOST_CHECK_EQUAL(orderedBarcodes[2][1].size(), 1);
  BOOST_CHECK(orderedBarcodes[2][1][0] == (std::pair<T, T>(5, 6)));
  BOOST_CHECK_EQUAL(orderedBarcodes[2][2].size(), 0);

  barcodes = persistence_on_slices<T>(s, {{0, 1, 1}, {2, 4, 2}, {2, 1, 0}}, {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}, false);
  std::vector<Test_multi_dimensional_barcode> orderedBarcodes2(barcodes.size());
  i = 0;
  for (const auto& dgm : barcodes) {
    orderedBarcodes2[i] = get_barcode(dgm);
    ++i;
  }

  BOOST_CHECK(orderedBarcodes == orderedBarcodes2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_batch_persistence, Slicer, list_of_tested_slicer_variants) {
  using Fil = typename Slicer::Filtration_value;

  auto cpx = build_simple_input_complex<Fil>();

  Slicer s(cpx);
  test_slicer_batch_persistence(s);
  Thread_safe_slicer wc = s.weak_copy();
  test_slicer_batch_persistence(wc);
  Thread_safe_slicer tss(s);
  test_slicer_batch_persistence(tss);
}

template <class Fil>
Multi_parameter_filtered_complex<Fil, I, D> build_simple_input_complex2() {
  using Complex = Multi_parameter_filtered_complex<Fil, I, D>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<T>;

  BC bc = {{},     {},     {},     {},     {},     {},     {},     {},     {},     {},           {0, 2},      {0, 3},
           {0, 5}, {0, 6}, {0, 8}, {1, 3}, {1, 9}, {2, 3}, {3, 7}, {4, 5}, {5, 6}, {10, 11, 17}, {12, 13, 20}};
  DC dc = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2};
  FC fc = {ini{0, 3},  ini{0, 2}, ini{0, 5}, ini{0, 2}, ini{0, 1},  ini{0, 0}, ini{0, 6}, ini{0, 6},
           ini{0, 4},  ini{0, 6}, ini{4, 5}, ini{9, 3}, ini{11, 3}, ini{8, 6}, ini{3, 4}, ini{2, 2},
           ini{10, 6}, ini{5, 5}, ini{7, 6}, ini{1, 1}, ini{6, 6},  ini{9, 5}, ini{11, 6}};

  return Complex(bc, dc, fc);
}

template <class Slicer>
void test_slicer_landscapes(Slicer& s) {
  std::vector<T> xGrid = {0, 2, 4, 6, 8, 10, 12};
  std::vector<T> yGrid = {0, 3, 6, 9, 12, 15};
  std::vector<T> direction = {2, 1};
  std::vector<T> ks = {0, 1, 2};

  std::vector<double> landscapes = compute_slicer_landscapes_on_grid(s, xGrid, yGrid, direction, 2, 3, 1, 0, ks, true);
  BOOST_CHECK_EQUAL(landscapes.size(), ks.size() * xGrid.size() * yGrid.size());
  Gudhi::Simple_mdspan view(landscapes.data(), ks.size(), xGrid.size(), yGrid.size());

  std::vector<std::vector<std::vector<double>>> realLS = {{{0., 0., 0., 0., 0., 0.},
                                                           {0., 1., 1., 1., 1., 1.},
                                                           {0., 2., 2., 1., 1., 1.},
                                                           {0., 3., 3., 1., 2., 2.},
                                                           {0., 3., 4., 1., 3., 3.},
                                                           {0., 3., 5., 1., 4., 4.},
                                                           {0., 3., 6., 1., 4., 5.}},

                                                          {{0., 0., 0., 0., 0., 0.},
                                                           {0., 1., 1., 1., 1., 1.},
                                                           {0., 1., 2., 0., 1., 1.},
                                                           {0., 1., 1., 0., 2., 2.},
                                                           {0., 1., 0., 0., 2., 1.},
                                                           {0., 0.5, 0., 0., 1.5, 1.},
                                                           {0., 0., 0., 0., 0.5, 0.}},

                                                          {{0., 0., 0., 0., 0., 0.},
                                                           {0., 0., 1., 1., 1., 1.},
                                                           {0., 0., 0.5, 0., 1., 1.},
                                                           {0., 0., 0., 0., 1., 1.},
                                                           {0., 0., 0., 0., 1., 1.},
                                                           {0., 0., 0., 0., 0.5, 0.},
                                                           {0., 0., 0., 0., 0., 0.}}};

  for (std::size_t i = 0; i < ks.size(); ++i) {
    for (std::size_t j = 0; j < xGrid.size(); ++j) {
      for (std::size_t k = 0; k < yGrid.size(); ++k) {
        BOOST_CHECK_EQUAL(view(i, j, k), realLS[i][j][k]);
      }
    }
  }
}
