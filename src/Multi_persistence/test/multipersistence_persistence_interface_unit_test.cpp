/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024-25 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstdint>  //std::uint32_t

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/Persistence_interface_matrix.h>
#include <gudhi/Multi_persistence/Persistence_interface_cohomology.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_parameter_filtration.h>

using Gudhi::multi_filtration::Multi_parameter_filtration;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;
using Gudhi::multi_persistence::Persistence_interface_matrix;
using Gudhi::multi_persistence::Persistence_interface_cohomology;

using Fil = Multi_parameter_filtration<double>;
using Complex = Multi_parameter_filtered_complex<Fil>;
using FC = typename Complex::Filtration_value_container;
using BC = typename Complex::Boundary_container;
using DC = typename Complex::Dimension_container;

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

typedef boost::mpl::list<Persistence_interface_matrix<Multi_persistence_r_options>,
                         Persistence_interface_matrix<Multi_persistence_ru_options>,
                         Persistence_interface_matrix<Multi_persistence_chain_options>,
                         Persistence_interface_cohomology<Fil>>
    list_of_tested_variants;

Complex build_complex()
{
  using ini = std::initializer_list<typename Fil::value_type>;

  BC bc = {{}, {0, 5}, {}, {}, {}, {}, {0, 4}, {5, 3}, {5, 4}, {2, 3}, {2, 4}, {3, 4}, {1, 6, 8}, {9, 10, 11}};
  DC dc = {0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2};
  FC fc = {ini{0, 1, 2},
           ini{0, 1, 2},
           ini{0, 1, 2},
           ini{3, 4, 5},
           ini{3, 4, 5},
           ini{6, 7, 8},
           ini{6, 7, 8},
           ini{6, 7, 8},
           ini{6, 7, 8},
           ini{6, 7, 8},
           ini{6, 7, 8},
           ini{6, 7, 8},
           ini{6, 7, 8},
           ini{6, 7, 8}};

  return Complex(bc, dc, fc);
}

Complex build_complex_vine()
{
  using ini = std::initializer_list<typename Fil::value_type>;

  BC bc = {{}, {}, {}, {0, 1}, {1, 2}, {0, 2}, {3, 4, 5}, {}, {1, 7}};
  DC dc = {0, 0, 0, 1, 1, 1, 2, 0, 1};
  FC fc = {ini{0, 1, 2},
           ini{0, 1, 2},
           ini{0, 1, 2},
           ini{3, 4, 5},
           ini{3, 4, 5},
           ini{6, 7, 8},
           ini{6, 7, 8},
           ini{6, 7, 8},
           ini{6, 7, 8}};

  return Complex(bc, dc, fc);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(interface_constructors, Interface, list_of_tested_variants)
{
  using Map = typename Interface::Map;

  Interface empty;
  BOOST_CHECK(!empty.is_initialized());

  Complex cpx = build_complex();
  Map permutation = {0, 5, 2, 3, 4, 6, 8, 1, 10, 11, 9, 7, 12, 13};

  Interface inter(cpx, permutation);
  BOOST_CHECK(inter.is_initialized());

  for (int i = 0; i < 14; ++i) {
    if (i < 5)
      BOOST_CHECK_EQUAL(inter.get_dimension(i), 0);
    else if (i < 12)
      BOOST_CHECK_EQUAL(inter.get_dimension(i), 1);
    else
      BOOST_CHECK_EQUAL(inter.get_dimension(i), 2);
  }

  Interface copy(inter, permutation);
  BOOST_CHECK(copy.is_initialized());

  for (int i = 0; i < 14; ++i) {
    if (i < 5)
      BOOST_CHECK_EQUAL(copy.get_dimension(i), 0);
    else if (i < 12)
      BOOST_CHECK_EQUAL(copy.get_dimension(i), 1);
    else
      BOOST_CHECK_EQUAL(copy.get_dimension(i), 2);
  }

  Interface move(std::move(inter), permutation);
  BOOST_CHECK(move.is_initialized());
  BOOST_CHECK(!inter.is_initialized());

  for (int i = 0; i < 14; ++i) {
    if (i < 5)
      BOOST_CHECK_EQUAL(move.get_dimension(i), 0);
    else if (i < 12)
      BOOST_CHECK_EQUAL(move.get_dimension(i), 1);
    else
      BOOST_CHECK_EQUAL(move.get_dimension(i), 2);
  }
}

template <class Barcode, class Bar>
void test_barcode_equality(const Barcode& barcode1, const std::vector<Bar>& barcode2)
{
  BOOST_CHECK_EQUAL(barcode1.size(), barcode2.size());
  std::vector<Bar> sorted;
  for (const Bar& b : barcode1) sorted.push_back(b);
  std::sort(sorted.begin(), sorted.end(), [](const Bar& b1, const Bar& b2){
    return b1.birth < b2.birth;
  });
  for (unsigned int i = 0; i < barcode2.size(); ++i) {
    const auto& b1 = sorted[i];
    const auto& b2 = barcode2[i];
    BOOST_CHECK_EQUAL(b1.dim, b2.dim);
    BOOST_CHECK_EQUAL(b1.birth, b2.birth);
    BOOST_CHECK_EQUAL(b1.death, b2.death);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(interface_barcode, Interface, list_of_tested_variants)
{
  using Map = typename Interface::Map;
  using Bar = typename Interface::Bar;

  Complex cpx = build_complex();
  Map permutation = {0, 5, 2, 3, 4, 6, 8, 1, 10, 11, 9, 7, 12, 13};
  std::vector<Bar> realBarcode;

  realBarcode = {Bar(0, Bar::inf, 0),
                 Bar(1, 12, 1),
                 Bar(2, 10, 0),
                 Bar(3, 11, 0),
                 Bar(4, 6, 0),
                 Bar(5, 8, 0),
                 Bar(7, Bar::inf, 1),
                 Bar(9, 13, 1)};

  Interface inter(cpx, permutation);
  auto barcode = inter.get_barcode();

  test_barcode_equality(barcode, realBarcode);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(interface_vine, Interface, list_of_tested_variants)
{
  using Map = typename Interface::Map;
  using Bar = typename Interface::Bar;

  if constexpr (Interface::is_vine) {
    Complex cpx = build_complex_vine();
    Map permutation = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Bar> realBarcode = {Bar(0, Bar::inf, 0), Bar(1, 3, 0), Bar(2, 4, 0), Bar(5, 6, 1), Bar(7, 8, 0)};

    Interface inter(cpx, permutation);
    auto barcode = inter.get_barcode();
    test_barcode_equality(barcode, realBarcode);

    inter.vine_swap(6);
    std::swap(permutation[6], permutation[7]);
    inter.vine_swap(5);
    std::swap(permutation[5], permutation[6]);
    inter.vine_swap(4);
    std::swap(permutation[4], permutation[5]);
    inter.vine_swap(3);
    std::swap(permutation[3], permutation[4]);
    inter.vine_swap(7);
    std::swap(permutation[7], permutation[8]);
    barcode = inter.get_barcode();
    test_barcode_equality(barcode, realBarcode);

    inter.vine_swap(0);
    std::swap(permutation[0], permutation[1]);
    realBarcode[0].birth = 1;
    realBarcode[1].birth = 0;
    std::swap(realBarcode[0], realBarcode[1]);  // reorder by birth
    barcode = inter.get_barcode();
    test_barcode_equality(barcode, realBarcode);

    inter.vine_swap(4);
    std::swap(permutation[4], permutation[5]);
    barcode = inter.get_barcode();
    test_barcode_equality(barcode, realBarcode);

    inter.vine_swap(5);
    std::swap(permutation[5], permutation[6]);
    realBarcode[0].death = 5;
    realBarcode[3].birth = 3;
    barcode = inter.get_barcode();
    test_barcode_equality(barcode, realBarcode);

    inter.vine_swap(6);
    std::swap(permutation[6], permutation[7]);
    barcode = inter.get_barcode();
    test_barcode_equality(barcode, realBarcode);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(interface_rep_cycles, Interface, list_of_tested_variants)
{
  using Map = typename Interface::Map;

  if constexpr (Interface::has_rep_cycles) {
    using Cy = typename Interface::Cycle;

    Complex cpx = build_complex();
    Map permutation = {0, 5, 2, 3, 4, 6, 8, 1, 10, 11, 9, 7, 12, 13};

    Interface inter(cpx, permutation);
    const auto& cycles = inter.get_representative_cycles(true);
    BOOST_CHECK_EQUAL(cycles.size(), 8);

    BOOST_CHECK((cycles(0) == Cy{0}));
    BOOST_CHECK((cycles(1) == Cy{0, 5}));
    BOOST_CHECK((cycles(2) == Cy{0, 2}));
    BOOST_CHECK((cycles(3) == Cy{0, 3}));
    BOOST_CHECK((cycles(4) == Cy{0, 4}));
    BOOST_CHECK((cycles(5) == Cy{6, 8, 1}));
    BOOST_CHECK((cycles(6) == Cy{10, 11, 9}));
    BOOST_CHECK((cycles(7) == Cy{8, 11, 7}));

    BOOST_CHECK((inter.get_representative_cycle(0, false) == Cy{0}));
    BOOST_CHECK((inter.get_representative_cycle(1, false) == Cy{0, 5}));
    BOOST_CHECK((inter.get_representative_cycle(2, false) == Cy{0, 2}));
    BOOST_CHECK((inter.get_representative_cycle(3, false) == Cy{0, 3}));
    BOOST_CHECK((inter.get_representative_cycle(4, false) == Cy{0, 4}));
    BOOST_CHECK((inter.get_representative_cycle(5, false) == Cy{6, 8, 1}));
    BOOST_CHECK((inter.get_representative_cycle(6, false) == Cy{10, 11, 9}));
    BOOST_CHECK((inter.get_representative_cycle(7, false) == Cy{8, 11, 7}));
  }
}
