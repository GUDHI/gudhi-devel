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

#include <gudhi/vineyard_builder.h>

#include "vy_test_utilities.h"

using Gudhi::vineyard::Vineyard_builder;

template <typename T>
using Vine = Gudhi::vineyard::Vine<T, int>;
using Cycle = std::vector<int>;

using option_list = boost::mpl::list<Chain_vineyard_options, RU_vineyard_options>;

template <class VB>
std::vector<Vine<typename VB::value_type> > get_vineyard(const VB& vb)
{
  using V = Vine<typename VB::value_type>;

  std::vector<V> res;

  if constexpr (VB::has_flat_vineyard()) {
    const auto& numVinesByDim = vb.get_number_of_vines_by_dimension();
    unsigned int i = 0;

    const auto& vy = vb.get_current_vineyard();
    for (unsigned int d = 0; d < vy.size(); ++d) {
      auto nv = numVinesByDim[d];
      auto numUpdates = vy[d].size() / nv;
      res.resize(i + nv, V(d, numUpdates));
      for (unsigned int k = 0; k < numUpdates; ++k) {
        for (unsigned int l = 0; l < nv; ++l) {
          const auto& p = vy[d][l + k * nv];
          res[i + l].add_pair(p[0], p[1]);
        }
      }
      i = res.size();
    }
  } else {
    res = vb.get_current_vineyard();
  }

  std::sort(res.begin(), res.end(), [](const V& v1, const V& v2) {
    if (v1.size() == 0) return false;  // both have the same size
    if (v1.get_dimension() != v2.get_dimension()) return v1.get_dimension() < v2.get_dimension();
    return v1.get_pairs() < v2.get_pairs();
  });

  return res;
}

template <class VB>
std::vector<std::tuple<Cycle, typename VB::Dimension, typename VB::Index> > get_ordered_cycles(const VB& vb)
{
  using C = std::tuple<Cycle, typename VB::Dimension, typename VB::Index>;
  const auto& cs = vb.get_latest_representative_cycles();
  std::vector<C> res(cs.size());
  for (unsigned int i = 0; i < cs.size(); ++i) {
    res[i] = C(get_cycle(get<0>(cs[i])), get<1>(cs[i]), get<2>(cs[i]));
  }
  std::sort(res.begin(), res.end(), [](const C& c1, const C& c2) {
    if (get<1>(c1) != get<1>(c2)) return get<1>(c1) < get<1>(c2);
    return get<0>(c1) < get<0>(c2);
  });
  return res;
}

template <class VB>
void test_initialization()
{
  using T = typename VB::value_type;
  using Bar = typename VB::Bar;
  using P = typename Vine<T>::Coordinate;

  auto [bc, dc, fc] = build_simple_input_complex<T>();

  VB vb_no_c;
  VB vb_c(true);
  VB vb_c_dim1(true, 1);

  auto vy = get_vineyard(vb_no_c);
  BOOST_CHECK(vy.empty());
  vy = get_vineyard(vb_c);
  BOOST_CHECK(vy.empty());
  vy = get_vineyard(vb_c_dim1);
  BOOST_CHECK(vy.empty());

  vb_no_c.initialize(bc, dc, fc);
  vb_c.initialize(bc, dc, fc);
  vb_c_dim1.initialize(bc, dc, fc);

  vy = get_vineyard(vb_no_c);
  BOOST_CHECK_EQUAL(vy.size(), 5);
  BOOST_CHECK_EQUAL(vy[0].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[0].size(), 1);
  BOOST_CHECK(vy[0].get_pair(0) == (P{1, 3}));
  BOOST_CHECK_EQUAL(vy[1].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[1].size(), 1);
  BOOST_CHECK(vy[1].get_pair(0) == (P{1, Bar::inf}));
  BOOST_CHECK_EQUAL(vy[2].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[2].size(), 1);
  BOOST_CHECK(vy[2].get_pair(0) == (P{2, 4}));
  BOOST_CHECK_EQUAL(vy[3].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[3].size(), 1);
  BOOST_CHECK(vy[3].get_pair(0) == (P{Bar::inf, Bar::inf}));
  BOOST_CHECK_EQUAL(vy[4].get_dimension(), 1);
  BOOST_CHECK_EQUAL(vy[4].size(), 1);
  BOOST_CHECK(vy[4].get_pair(0) == (P{6, Bar::inf}));

  BOOST_CHECK(vy == get_vineyard(vb_c));
  BOOST_CHECK(vy == get_vineyard(vb_c_dim1));

  BOOST_CHECK_THROW(vb_no_c.get_latest_representative_cycles(), std::invalid_argument);

  auto cycles = get_ordered_cycles(vb_c);
  if constexpr (VB::has_flat_vineyard()) {
    BOOST_CHECK_EQUAL(get<2>(cycles[3]), 0);
  } else {
    BOOST_CHECK_EQUAL(get<2>(cycles[3]), 4);
  }

  cycles = get_ordered_cycles(vb_c_dim1);
  if constexpr (VB::has_flat_vineyard()) {
    BOOST_CHECK_EQUAL(get<2>(cycles[0]), 0);
  } else {
    BOOST_CHECK_EQUAL(get<2>(cycles[0]), 4);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(vyb_initialization, Option, option_list)
{
  test_initialization<Vineyard_builder<double, Option, false> >();
  test_initialization<Vineyard_builder<double, Option, true> >();
  test_initialization<Vineyard_builder<int, Option, false> >();
  test_initialization<Vineyard_builder<int, Option, true> >();
  test_initialization<Vineyard_builder<unsigned int, Option, false> >();
  test_initialization<Vineyard_builder<unsigned int, Option, true> >();
}

template <class VB>
void test_update(){
  using T = typename VB::value_type;
  using Bar = typename VB::Bar;
  using P = typename Vine<T>::Coordinate;

  auto [bc, dc, fc] = build_simple_input_complex<T>();

  VB vb_no_c;
  VB vb_c(true);
  VB vb_c_dim1(true, 1);

  auto vy = get_vineyard(vb_no_c);
  BOOST_CHECK(vy.empty());
  vy = get_vineyard(vb_c);
  BOOST_CHECK(vy.empty());
  vy = get_vineyard(vb_c_dim1);
  BOOST_CHECK(vy.empty());

  vb_no_c.initialize(bc, dc, fc);
  vb_c.initialize(bc, dc, fc);
  vb_c_dim1.initialize(bc, dc, fc);

  vb_no_c.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});
  vb_c.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});
  vb_c_dim1.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});

  vy = get_vineyard(vb_no_c);
  BOOST_CHECK_EQUAL(vy.size(), 5);
  BOOST_CHECK_EQUAL(vy[0].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[0].size(), 2);
  BOOST_CHECK(vy[0].get_pair(0) == (P{1, 3}));
  BOOST_CHECK(vy[0].get_pair(1) == (P{2, 4}));
  BOOST_CHECK_EQUAL(vy[1].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[1].size(), 2);
  BOOST_CHECK(vy[1].get_pair(0) == (P{1, Bar::inf}));
  BOOST_CHECK(vy[1].get_pair(1) == (P{1, Bar::inf}));
  BOOST_CHECK_EQUAL(vy[2].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[2].size(), 2);
  BOOST_CHECK(vy[2].get_pair(0) == (P{2, 4}));
  BOOST_CHECK(vy[2].get_pair(1) == (P{1, 3}));
  BOOST_CHECK_EQUAL(vy[3].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[3].size(), 2);
  BOOST_CHECK(vy[3].get_pair(0) == (P{Bar::inf, Bar::inf}));
  BOOST_CHECK(vy[3].get_pair(1) == (P{7, 7}));
  BOOST_CHECK_EQUAL(vy[4].get_dimension(), 1);
  BOOST_CHECK_EQUAL(vy[4].size(), 2);
  BOOST_CHECK(vy[4].get_pair(0) == (P{6, Bar::inf}));
  BOOST_CHECK(vy[4].get_pair(1) == (P{6, 6}));

  BOOST_CHECK(vy == get_vineyard(vb_c));
  BOOST_CHECK(vy == get_vineyard(vb_c_dim1));

  vb_no_c.update(FC<double>{0, -1, 1, 1, 3, 4, 4, 6, 6});
  vb_c.update(FC<double>{0, -1, 1, 1, 3, 4, 4, 6, 6});
  vb_c_dim1.update(FC<double>{0, -1, 1, 1, 3, 4, 4, 6, 6});

  vy = get_vineyard(vb_no_c);
  BOOST_CHECK_EQUAL(vy.size(), 5);
  BOOST_CHECK_EQUAL(vy[0].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[0].size(), 3);
  BOOST_CHECK(vy[0].get_pair(0) == (P{1, 3}));
  BOOST_CHECK(vy[0].get_pair(1) == (P{2, 4}));
  BOOST_CHECK(vy[0].get_pair(2) == (P{1, 3}));
  BOOST_CHECK_EQUAL(vy[1].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[1].size(), 3);
  BOOST_CHECK(vy[1].get_pair(0) == (P{1, Bar::inf}));
  BOOST_CHECK(vy[1].get_pair(1) == (P{1, Bar::inf}));
  BOOST_CHECK(vy[1].get_pair(2) == (P{-1, Bar::inf}));
  BOOST_CHECK_EQUAL(vy[2].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[2].size(), 3);
  BOOST_CHECK(vy[2].get_pair(0) == (P{2, 4}));
  BOOST_CHECK(vy[2].get_pair(1) == (P{1, 3}));
  BOOST_CHECK(vy[2].get_pair(2) == (P{0, 1}));
  BOOST_CHECK_EQUAL(vy[3].get_dimension(), 0);
  BOOST_CHECK_EQUAL(vy[3].size(), 3);
  BOOST_CHECK(vy[3].get_pair(0) == (P{Bar::inf, Bar::inf}));
  BOOST_CHECK(vy[3].get_pair(1) == (P{7, 7}));
  BOOST_CHECK(vy[3].get_pair(2) == (P{6, 6}));
  BOOST_CHECK_EQUAL(vy[4].get_dimension(), 1);
  BOOST_CHECK_EQUAL(vy[4].size(), 3);
  BOOST_CHECK(vy[4].get_pair(0) == (P{6, Bar::inf}));
  BOOST_CHECK(vy[4].get_pair(1) == (P{6, 6}));
  BOOST_CHECK(vy[4].get_pair(2) == (P{4, 4}));

  BOOST_CHECK(vy == get_vineyard(vb_c));
  BOOST_CHECK(vy == get_vineyard(vb_c_dim1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(vyb_update, Option, option_list)
{
  test_update<Vineyard_builder<double, Option, false> >();
  test_update<Vineyard_builder<double, Option, true> >();
  test_update<Vineyard_builder<int, Option, false> >();
  test_update<Vineyard_builder<int, Option, true> >();
}

template<class Option>
void test_rep_cycles(){}

template <>
void test_rep_cycles<Chain_vineyard_options>()
{
  using VB = Vineyard_builder<double, Chain_vineyard_options, false>;
  using T = typename VB::value_type;

  auto [bc, dc, fc] = build_simple_input_complex<T>();

  VB vb_c(true);
  VB vb_c_dim1(true, 1);

  vb_c.initialize(bc, dc, fc);
  vb_c_dim1.initialize(bc, dc, fc);

  // TODO: instead of testing the value of get<2>(...), we should test if that value corresponds
  // to the right bar in vb.get_current_vineyard(), as the indices can potentially change with modification of
  // the algorithms without being wrong. But that is also true for other things tested here and is quite bothersome.
  // For now, just be careful about it.
  auto cycles = get_ordered_cycles(vb_c);
  BOOST_CHECK_EQUAL(cycles.size(), 4);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{0}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 0);
  BOOST_CHECK(get<0>(cycles[1]) == (Cycle{0, 1}));
  BOOST_CHECK_EQUAL(get<1>(cycles[1]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[1]), 2);
  BOOST_CHECK(get<0>(cycles[2]) == (Cycle{0, 2}));
  BOOST_CHECK_EQUAL(get<1>(cycles[2]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[2]), 1);
  BOOST_CHECK(get<0>(cycles[3]) == (Cycle{3, 4, 5}));
  BOOST_CHECK_EQUAL(get<1>(cycles[3]), 1);
  BOOST_CHECK_EQUAL(get<2>(cycles[3]), 4);

  cycles = get_ordered_cycles(vb_c_dim1);
  BOOST_CHECK_EQUAL(cycles.size(), 1);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{3, 4, 5}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 1);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 4);

  vb_c.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});
  vb_c_dim1.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});

  cycles = get_ordered_cycles(vb_c);
  BOOST_CHECK_EQUAL(cycles.size(), 3);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{0}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 0);
  BOOST_CHECK(get<0>(cycles[1]) == (Cycle{0, 1}));
  BOOST_CHECK_EQUAL(get<1>(cycles[1]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[1]), 2);
  BOOST_CHECK(get<0>(cycles[2]) == (Cycle{1, 2}));
  BOOST_CHECK_EQUAL(get<1>(cycles[2]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[2]), 1);

  cycles = get_ordered_cycles(vb_c_dim1);
  BOOST_CHECK_EQUAL(cycles.size(), 0);

  vb_c.update(FC<double>{0, -1, 1, 1, 3, 4, 5, 6, 7});
  vb_c_dim1.update(FC<double>{0, -1, 1, 1, 3, 4, 5, 6, 7});
  
  cycles = get_ordered_cycles(vb_c);
  BOOST_CHECK_EQUAL(cycles.size(), 5);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{0, 1}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 2);
  BOOST_CHECK(get<0>(cycles[1]) == (Cycle{0, 7}));
  BOOST_CHECK_EQUAL(get<1>(cycles[1]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[1]), 3);
  BOOST_CHECK(get<0>(cycles[2]) == (Cycle{1}));
  BOOST_CHECK_EQUAL(get<1>(cycles[2]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[2]), 0);
  BOOST_CHECK(get<0>(cycles[3]) == (Cycle{1, 2}));
  BOOST_CHECK_EQUAL(get<1>(cycles[3]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[3]), 1);
  BOOST_CHECK(get<0>(cycles[4]) == (Cycle{3, 4, 5}));
  BOOST_CHECK_EQUAL(get<1>(cycles[4]), 1);
  BOOST_CHECK_EQUAL(get<2>(cycles[4]), 4);

  cycles = get_ordered_cycles(vb_c_dim1);
  BOOST_CHECK_EQUAL(cycles.size(), 1);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{3, 4, 5}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 1);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 4);
}

template <>
void test_rep_cycles<RU_vineyard_options>()
{
  using VB = Vineyard_builder<double, RU_vineyard_options, false>;
  using T = typename VB::value_type;

  auto [bc, dc, fc] = build_simple_input_complex<T>();

  VB vb_c(true);
  VB vb_c_dim1(true, 1);

  vb_c.initialize(bc, dc, fc);
  vb_c_dim1.initialize(bc, dc, fc);

  auto cycles = get_ordered_cycles(vb_c);
  BOOST_CHECK_EQUAL(cycles.size(), 4);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{0}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 0);
  BOOST_CHECK(get<0>(cycles[1]) == (Cycle{0, 2}));
  BOOST_CHECK_EQUAL(get<1>(cycles[1]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[1]), 1);
  BOOST_CHECK(get<0>(cycles[2]) == (Cycle{1, 2}));
  BOOST_CHECK_EQUAL(get<1>(cycles[2]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[2]), 2);
  BOOST_CHECK(get<0>(cycles[3]) == (Cycle{3, 4, 5}));
  BOOST_CHECK_EQUAL(get<1>(cycles[3]), 1);
  BOOST_CHECK_EQUAL(get<2>(cycles[3]), 4);

  cycles = get_ordered_cycles(vb_c_dim1);
  BOOST_CHECK_EQUAL(cycles.size(), 1);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{3, 4, 5}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 1);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 4);

  vb_c.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});
  vb_c_dim1.update(FC<double>{1, 1, 2, 3, 4, 6, 6, 7, 7});

  cycles = get_ordered_cycles(vb_c);
  BOOST_CHECK_EQUAL(cycles.size(), 3);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{0}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 0);
  BOOST_CHECK(get<0>(cycles[1]) == (Cycle{0, 1}));
  BOOST_CHECK_EQUAL(get<1>(cycles[1]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[1]), 2);
  BOOST_CHECK(get<0>(cycles[2]) == (Cycle{1, 2}));
  BOOST_CHECK_EQUAL(get<1>(cycles[2]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[2]), 1);

  cycles = get_ordered_cycles(vb_c_dim1);
  BOOST_CHECK_EQUAL(cycles.size(), 0);

  vb_c.update(FC<double>{0, -1, 1, 1, 3, 4, 5, 6, 7});
  vb_c_dim1.update(FC<double>{0, -1, 1, 1, 3, 4, 5, 6, 7});
  
  cycles = get_ordered_cycles(vb_c);
  BOOST_CHECK_EQUAL(cycles.size(), 5);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{0, 1}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 2);
  BOOST_CHECK(get<0>(cycles[1]) == (Cycle{1}));
  BOOST_CHECK_EQUAL(get<1>(cycles[1]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[1]), 0);
  BOOST_CHECK(get<0>(cycles[2]) == (Cycle{1, 2}));
  BOOST_CHECK_EQUAL(get<1>(cycles[2]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[2]), 1);
  BOOST_CHECK(get<0>(cycles[3]) == (Cycle{1, 7}));
  BOOST_CHECK_EQUAL(get<1>(cycles[3]), 0);
  BOOST_CHECK_EQUAL(get<2>(cycles[3]), 3);
  BOOST_CHECK(get<0>(cycles[4]) == (Cycle{3, 4, 5}));
  BOOST_CHECK_EQUAL(get<1>(cycles[4]), 1);
  BOOST_CHECK_EQUAL(get<2>(cycles[4]), 4);

  cycles = get_ordered_cycles(vb_c_dim1);
  BOOST_CHECK_EQUAL(cycles.size(), 1);
  BOOST_CHECK(get<0>(cycles[0]) == (Cycle{3, 4, 5}));
  BOOST_CHECK_EQUAL(get<1>(cycles[0]), 1);
  BOOST_CHECK_EQUAL(get<2>(cycles[0]), 4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(vy_rep_cycles, Option, option_list) {
  test_rep_cycles<Option>();
}
