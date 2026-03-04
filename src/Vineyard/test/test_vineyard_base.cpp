/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "vineyard"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/vineyard_base.h>

#include "vy_test_utilities.h"

using option_list = boost::mpl::list<Chain_vineyard_options, RU_vineyard_options>;

template <class Bar, class Barcode>
std::vector<Bar> get_barcode(const Barcode& bc)
{
  std::vector<Bar> barcode(bc.begin(), bc.end());
  std::sort(barcode.begin(), barcode.end(), [](const Bar& b1, const Bar& b2) {
    if (b1.dim == b2.dim) return b1.birth < b2.birth;
    return b1.dim < b2.dim;
  });
  return barcode;
}

template <class Cycles>
std::vector<std::vector<int>> get_all_cycles(const Cycles& cs)
{
  std::vector<std::vector<int>> cycles(cs.size());
  unsigned int i = 0;
  for (const auto& c : cs) {
    cycles[i] = get_cycle(c);
    ++i;
  }
  std::sort(cycles.begin(), cycles.end());
  return cycles;
}

template <class V>
std::vector<std::vector<int>> get_all_cycles_individually(V& vy)
{
  std::vector<std::vector<int>> cycles(vy.get_current_barcode().size());
  for (unsigned int i = 0; i < cycles.size(); ++i) {
    cycles[i] = get_cycle(vy.get_current_representative_cycle(i));
  }
  std::sort(cycles.begin(), cycles.end());
  return cycles;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(vy_initialization, Option, option_list) {
  using V = Gudhi::vineyard::Vineyard_base<Option>;
  using Bar = typename V::Bar;

  auto [bc, dc, fc] = build_simple_input_complex();

  V empty;
  BOOST_CHECK(!empty.is_initialized());
  BOOST_CHECK_EQUAL(empty.get_current_barcode().size(), 0);

  V vy(bc, dc, fc);
  BOOST_CHECK(vy.is_initialized());

  auto order = vy.get_current_order();
  BOOST_CHECK_EQUAL(order.size(), 9);
  BOOST_CHECK_EQUAL(order[0], 0);
  BOOST_CHECK_EQUAL(order[1], 2);
  BOOST_CHECK_EQUAL(order[2], 1);
  BOOST_CHECK_EQUAL(order[3], 7);
  BOOST_CHECK_EQUAL(order[4], 5);
  BOOST_CHECK_EQUAL(order[5], 4);
  BOOST_CHECK_EQUAL(order[6], 3);
  BOOST_CHECK_EQUAL(order[7], 8);
  BOOST_CHECK_EQUAL(order[8], 6);

  auto barcode = get_barcode<Bar>(vy.get_current_barcode());
  BOOST_CHECK_EQUAL(barcode.size(), 5);
  BOOST_CHECK_EQUAL(barcode[0], Bar(0, Bar::inf, 0));
  BOOST_CHECK_EQUAL(barcode[1], Bar(1, 4, 0));
  BOOST_CHECK_EQUAL(barcode[2], Bar(2, 5, 0));
  BOOST_CHECK_EQUAL(barcode[3], Bar(7, 8, 0));
  BOOST_CHECK_EQUAL(barcode[4], Bar(3, 6, 1));

  auto cycles = get_all_cycles(vy.get_all_current_representative_cycles());
  BOOST_CHECK_EQUAL(cycles.size(), 5);

  empty.initialize(bc, dc, fc);
  BOOST_CHECK(order == empty.get_current_order());
  BOOST_CHECK(barcode == get_barcode<Bar>(empty.get_current_barcode()));
  BOOST_CHECK(cycles == get_all_cycles(empty.get_all_current_representative_cycles()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(vy_update, Option, option_list) {
  using V = Gudhi::vineyard::Vineyard_base<Option>;
  using Bar = typename V::Bar;

  auto [bc, dc, fc] = build_simple_input_complex();

  V vy(bc, dc, fc);

  vy.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});

  auto order = vy.get_current_order();
  BOOST_CHECK_EQUAL(order.size(), 9);
  BOOST_CHECK_EQUAL(order[0], 0);
  BOOST_CHECK_EQUAL(order[1], 1);
  BOOST_CHECK_EQUAL(order[2], 2);
  BOOST_CHECK_EQUAL(order[3], 7);
  BOOST_CHECK_EQUAL(order[4], 3);
  BOOST_CHECK_EQUAL(order[5], 4);
  BOOST_CHECK_EQUAL(order[6], 5);
  BOOST_CHECK_EQUAL(order[7], 8);
  BOOST_CHECK_EQUAL(order[8], 6);

  auto barcode = get_barcode<Bar>(vy.get_current_barcode());
  BOOST_CHECK_EQUAL(barcode.size(), 5);
  BOOST_CHECK_EQUAL(barcode[0], Bar(0, Bar::inf, 0));
  BOOST_CHECK_EQUAL(barcode[1], Bar(1, 3, 0));
  BOOST_CHECK_EQUAL(barcode[2], Bar(2, 4, 0));
  BOOST_CHECK_EQUAL(barcode[3], Bar(7, 8, 0));
  BOOST_CHECK_EQUAL(barcode[4], Bar(5, 6, 1));

  vy.update(FC<double>{0, -1, 1, 1, 3, 4, 4, 6, 6});

  order = vy.get_current_order();
  BOOST_CHECK_EQUAL(order.size(), 9);
  BOOST_CHECK_EQUAL(order[0], 1);
  BOOST_CHECK_EQUAL(order[1], 0);
  BOOST_CHECK_EQUAL(order[2], 2);
  BOOST_CHECK_EQUAL(order[3], 7);
  BOOST_CHECK_EQUAL(order[4], 3);
  BOOST_CHECK_EQUAL(order[5], 4);
  BOOST_CHECK_EQUAL(order[6], 5);
  BOOST_CHECK_EQUAL(order[7], 8);
  BOOST_CHECK_EQUAL(order[8], 6);

  barcode = get_barcode<Bar>(vy.get_current_barcode());
  BOOST_CHECK_EQUAL(barcode.size(), 5);
  BOOST_CHECK_EQUAL(barcode[0], Bar(0, 3, 0));
  BOOST_CHECK_EQUAL(barcode[1], Bar(1, Bar::inf, 0));
  BOOST_CHECK_EQUAL(barcode[2], Bar(2, 4, 0));
  BOOST_CHECK_EQUAL(barcode[3], Bar(7, 8, 0));
  BOOST_CHECK_EQUAL(barcode[4], Bar(5, 6, 1));
}

template<class Option>
void test_rep_cycles(){}

template<>
void test_rep_cycles<Chain_vineyard_options>(){
  using V = Gudhi::vineyard::Vineyard_base<Chain_vineyard_options>;
  using C = std::vector<int>;

  auto [bc, dc, fc] = build_simple_input_complex();

  V vy(bc, dc, fc);

  auto cycles1 = get_all_cycles(vy.get_all_current_representative_cycles());
  BOOST_CHECK_EQUAL(cycles1.size(), 5);
  BOOST_CHECK(cycles1[0] == C{0});
  BOOST_CHECK(cycles1[1] == (C{0, 1}));
  BOOST_CHECK(cycles1[2] == (C{0, 2}));
  BOOST_CHECK(cycles1[3] == (C{0, 7}));
  BOOST_CHECK(cycles1[4] == (C{3, 4, 5}));

  auto cycles2 = get_all_cycles_individually(vy);
  BOOST_CHECK(cycles1 == cycles2);

  auto cycles3 = get_all_cycles(vy.get_all_current_representative_cycles(true, 1));
  BOOST_CHECK_EQUAL(cycles3.size(), 1);
  BOOST_CHECK(cycles3[0] == (C{3, 4, 5}));

  vy.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});

  cycles1 = get_all_cycles(vy.get_all_current_representative_cycles());
  BOOST_CHECK_EQUAL(cycles1.size(), 5);
  BOOST_CHECK(cycles1[0] == C{0});
  BOOST_CHECK(cycles1[1] == (C{0, 1}));
  BOOST_CHECK(cycles1[2] == (C{0, 7}));
  BOOST_CHECK(cycles1[3] == (C{1, 2}));
  BOOST_CHECK(cycles1[4] == (C{3, 4, 5}));

  cycles2 = get_all_cycles_individually(vy);
  BOOST_CHECK(cycles1 == cycles2);

  cycles3 = get_all_cycles(vy.get_all_current_representative_cycles(true, 1));
  BOOST_CHECK_EQUAL(cycles3.size(), 1);
  BOOST_CHECK(cycles3[0] == (C{3, 4, 5}));

  vy.update(FC<double>{0, -1, 1, 1, 3, 4, 4, 6, 6});

  cycles1 = get_all_cycles(vy.get_all_current_representative_cycles());
  BOOST_CHECK_EQUAL(cycles1.size(), 5);
  BOOST_CHECK(cycles1[0] == (C{0, 1}));
  BOOST_CHECK(cycles1[1] == (C{0, 7}));
  BOOST_CHECK(cycles1[2] == C{1});
  BOOST_CHECK(cycles1[3] == (C{1, 2}));
  BOOST_CHECK(cycles1[4] == (C{3, 4, 5}));

  cycles2 = get_all_cycles_individually(vy);
  BOOST_CHECK(cycles1 == cycles2);

  cycles3 = get_all_cycles(vy.get_all_current_representative_cycles(true, 1));
  BOOST_CHECK_EQUAL(cycles3.size(), 1);
  BOOST_CHECK(cycles3[0] == (C{3, 4, 5}));
}

template<>
void test_rep_cycles<RU_vineyard_options>(){
  using V = Gudhi::vineyard::Vineyard_base<RU_vineyard_options>;
  using C = std::vector<int>;

  auto [bc, dc, fc] = build_simple_input_complex();

  V vy(bc, dc, fc);

  auto cycles1 = get_all_cycles(vy.get_all_current_representative_cycles());
  BOOST_CHECK_EQUAL(cycles1.size(), 5);
  BOOST_CHECK(cycles1[0] == C{0});
  BOOST_CHECK(cycles1[1] == (C{0, 2}));
  BOOST_CHECK(cycles1[2] == (C{1, 2}));
  BOOST_CHECK(cycles1[3] == (C{1, 7}));
  BOOST_CHECK(cycles1[4] == (C{3, 4, 5}));

  auto cycles2 = get_all_cycles_individually(vy);
  BOOST_CHECK(cycles1 == cycles2);

  auto cycles3 = get_all_cycles(vy.get_all_current_representative_cycles(true, 1));
  BOOST_CHECK_EQUAL(cycles3.size(), 1);
  BOOST_CHECK(cycles3[0] == (C{3, 4, 5}));

  vy.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});

  cycles1 = get_all_cycles(vy.get_all_current_representative_cycles());
  BOOST_CHECK_EQUAL(cycles1.size(), 5);
  BOOST_CHECK(cycles1[0] == C{0});
  BOOST_CHECK(cycles1[1] == (C{0, 1}));
  BOOST_CHECK(cycles1[2] == (C{1, 2}));
  BOOST_CHECK(cycles1[3] == (C{1, 7}));
  BOOST_CHECK(cycles1[4] == (C{3, 4, 5}));

  cycles2 = get_all_cycles_individually(vy);
  BOOST_CHECK(cycles1 == cycles2);

  cycles3 = get_all_cycles(vy.get_all_current_representative_cycles(true, 1));
  BOOST_CHECK_EQUAL(cycles3.size(), 1);
  BOOST_CHECK(cycles3[0] == (C{3, 4, 5}));

  vy.update(FC<double>{0, -1, 1, 1, 3, 4, 4, 6, 6});

  cycles1 = get_all_cycles(vy.get_all_current_representative_cycles());
  BOOST_CHECK_EQUAL(cycles1.size(), 5);
  BOOST_CHECK(cycles1[0] == (C{0, 1}));
  BOOST_CHECK(cycles1[1] == C{1});
  BOOST_CHECK(cycles1[2] == (C{1, 2}));
  BOOST_CHECK(cycles1[3] == (C{1, 7}));
  BOOST_CHECK(cycles1[4] == (C{3, 4, 5}));

  cycles2 = get_all_cycles_individually(vy);
  BOOST_CHECK(cycles1 == cycles2);

  cycles3 = get_all_cycles(vy.get_all_current_representative_cycles(true, 1));
  BOOST_CHECK_EQUAL(cycles3.size(), 1);
  BOOST_CHECK(cycles3[0] == (C{3, 4, 5}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(vy_rep_cycles, Option, option_list) {
  test_rep_cycles<Option>();
}
