/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <array>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/Module.h>
#include <gudhi/Multi_persistence/Box.h>

using Gudhi::multi_persistence::Box;
using Gudhi::multi_persistence::Module;

using list_of_tested_variants = boost::mpl::list<double, float, int, unsigned int>;

BOOST_AUTO_TEST_CASE_TEMPLATE(module_constructors, T, list_of_tested_variants) {
  Module<T> empty;
  BOOST_CHECK_EQUAL(empty.size(), 0);
  BOOST_CHECK(empty.begin() == empty.end());
  BOOST_CHECK_EQUAL(empty.get_number_of_parameters(), Module<T>::template get_null_value<int>());
  BOOST_CHECK_EQUAL(empty.get_max_dimension(), Module<T>::template get_null_value<typename Module<T>::Dimension>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(module_add, T, list_of_tested_variants) {
  using S_t = typename Module<T>::Summand_t;
  using S_f = Module<float>::Summand_t;

  const int numParam = 2;
  S_t s1({1, 3, 3, 1}, {3, 8, 7, 6, 6, 99}, numParam, 0);
  S_t s2({2, 3, 4, 1}, {7, 6}, numParam, 2);
  S_t s3({10, 30, 30, 10}, {30, 80, 70, 60, 60, 99}, numParam, 0);
  S_t s4({20, 30, 40, 10}, {70, 60}, numParam, 2);
  S_t s5(numParam);

  Module<T> m;

  m.add_summand(s1);
  m.add_summand(s2);

  BOOST_CHECK_EQUAL(m.size(), 2);
  BOOST_CHECK_EQUAL(m.get_number_of_parameters(), numParam);
  BOOST_CHECK_EQUAL(m.get_max_dimension(), 2);

  auto it = m.begin();
  BOOST_CHECK(*it == m.get_summand(0));
  BOOST_CHECK(*it == s1);
  ++it;
  BOOST_CHECK(*it == m.get_summand(1));
  BOOST_CHECK(*it == s2);
  ++it;
  BOOST_CHECK(it == m.end());

  m.resize(3, numParam);
  BOOST_CHECK_EQUAL(m.size(), 3);

  Module<T> m2;
  m2.add_summand(1, s3);
  m2.add_summand(0, s4, 4);
  m.merge(m2);
  BOOST_CHECK_EQUAL(m.size(), 5);
  BOOST_CHECK_EQUAL(m.get_number_of_parameters(), numParam);
  BOOST_CHECK_EQUAL(m.get_max_dimension(), 4);

  it = m.begin();
  BOOST_CHECK(*it == m.get_summand(0));
  BOOST_CHECK(*it == s1);
  ++it;
  BOOST_CHECK(*it == m.get_summand(1));
  BOOST_CHECK(*it == s2);
  ++it;
  BOOST_CHECK(*it == m.get_summand(2));
  BOOST_CHECK(*it == s5);
  ++it;
  s4.set_dimension(4);
  BOOST_CHECK(*it == m.get_summand(3));
  BOOST_CHECK(*it == s4);
  ++it;
  BOOST_CHECK(*it == m.get_summand(4));
  BOOST_CHECK(*it == s3);
  ++it;
  BOOST_CHECK(it == m.end());

  Module<float> m3;
  m3.merge(m, 0);
  BOOST_CHECK_EQUAL(m3.size(), 2);
  BOOST_CHECK_EQUAL(m3.get_number_of_parameters(), numParam);
  BOOST_CHECK_EQUAL(m3.get_max_dimension(), 0);

  S_f s1f({1, 3, 3, 1}, {3, 8, 7, 6, 6, 99}, numParam, 0);
  S_f s3f({10, 30, 30, 10}, {30, 80, 70, 60, 60, 99}, numParam, 0);

  auto it3 = m3.begin();
  BOOST_CHECK(*it3 == m3.get_summand(0));
  BOOST_CHECK(*it3 == s1f);
  ++it3;
  BOOST_CHECK(*it3 == m3.get_summand(1));
  BOOST_CHECK(*it3 == s3f);
  ++it3;
  BOOST_CHECK(it3 == m3.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(module_modifiers, T, list_of_tested_variants) {
  using S = typename Module<T>::Summand_t;
  using B = typename S::Births;
  using D = typename S::Deaths;
  T inf = Module<T>::T_inf;

  const int numParam = 2;
  S s1({1, 6, 2, 5}, {8, 9, 12, 8}, numParam, 0);
  S s2({4, 3}, {inf, inf}, numParam, 2);
  S s3({5, 2, 8, 1}, {13, 7, 14, 6}, numParam, 0);
  S s5(B::inf(numParam), D::inf(numParam), 0);
  S s6(B::minus_inf(numParam), D::minus_inf(numParam), 2);

  Module<T> m;

  m.add_summand(s1);
  m.add_summand(s2);
  m.add_summand(s5);
  m.add_summand(s3);
  m.add_summand(s6);
  m.add_summand(s5);

  BOOST_CHECK_EQUAL(m.size(), 6);
  BOOST_CHECK(!m.get_summand(0).get_upset().is_plus_inf() && !m.get_summand(0).get_upset().is_minus_inf());
  BOOST_CHECK(!m.get_summand(1).get_upset().is_plus_inf() && !m.get_summand(1).get_upset().is_minus_inf());
  BOOST_CHECK(m.get_summand(2).get_upset().is_plus_inf());
  BOOST_CHECK(!m.get_summand(3).get_upset().is_plus_inf() && !m.get_summand(3).get_upset().is_minus_inf());
  BOOST_CHECK(m.get_summand(4).get_upset().is_minus_inf());
  BOOST_CHECK(m.get_summand(5).get_upset().is_plus_inf());

  m.clean();
  BOOST_CHECK_EQUAL(m.size(), 3);
  BOOST_CHECK(m.get_summand(0) == s1);
  BOOST_CHECK(m.get_summand(1) == s2);
  BOOST_CHECK(m.get_summand(2) == s3);

  S s1i({1, 5}, {8, 9, 12, 8}, numParam, 0);
  S s2i({4, 3}, {inf, inf}, numParam, 2);
  S s3i({5, 2, 8, 1}, {14, 7}, numParam, 0);

  m.fill(2);
  BOOST_CHECK_EQUAL(m.size(), 3);
  BOOST_CHECK(m.get_summand(0) == s1i);
  BOOST_CHECK(m.get_summand(1) == s2i);
  BOOST_CHECK(m.get_summand(2) == s3i);

  S s1r({2, 15}, {16, 27, 24, 24}, numParam, 0);
  S s2r({8, 9}, {inf, inf}, numParam, 2);
  S s3r({10, 6, 16, 3}, {28, 21}, numParam, 0);

  m.rescale(std::vector<int>{2, 3});
  BOOST_CHECK_EQUAL(m.size(), 3);
  BOOST_CHECK(m.get_summand(0) == s1r);
  BOOST_CHECK(m.get_summand(1) == s2r);
  BOOST_CHECK(m.get_summand(2) == s3r);

  S s1t({0, 12}, {14, 24, 22, 21}, numParam, 0);
  S s2t({6, 6}, {inf, inf}, numParam, 2);
  S s3t({8, 3, 14, 0}, {26, 18}, numParam, 0);

  m.translate(std::vector<int>{-2, -3});
  BOOST_CHECK_EQUAL(m.size(), 3);
  BOOST_CHECK(m.get_summand(0) == s1t);
  BOOST_CHECK(m.get_summand(1) == s2t);
  BOOST_CHECK(m.get_summand(2) == s3t);

  std::vector<std::vector<T>> grid = {{13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                                      {18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};
  S s1g({13, 6}, {inf, inf, inf, inf}, numParam, 0);
  S s2g({7, 12}, {inf, inf}, numParam, 2);
  S s3g({5, 15, inf, 18}, {inf, 0}, numParam, 0);

  m.evaluate_in_grid(grid);
  BOOST_CHECK_EQUAL(m.size(), 3);
  BOOST_CHECK(m.get_summand(0) == s1g);
  BOOST_CHECK(m.get_summand(1) == s2g);
  BOOST_CHECK(m.get_summand(2) == s3g);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(module_other, T, list_of_tested_variants) {
  using S_t = typename Module<T>::Summand_t;

  const int numParam = 2;
  S_t s1({1, 6}, {3, 8}, numParam, 0);
  S_t s2({2, 3}, {6, 10}, numParam, 2);
  S_t s3({4, 6}, {8, 11}, numParam, 0);
  S_t s4({4, 8}, {6, 9}, numParam, 2);

  Module<T> m;

  m.add_summand(s1);
  m.add_summand(s2);
  m.add_summand(s3);
  m.add_summand(s4);

  BOOST_CHECK(m.compute_bounds() == Box<T>({1, 3}, {8, 11}));

  std::vector<std::vector<std::array<double, 2>>> barcode({{{-0.5, 0.5}, {1, 3}}, {}, {{0, 2}, {1, 2}}});

  BOOST_CHECK(barcode == m.get_barcode_from_line({{2, 7}, {2, 1}}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(module_serialization, T, list_of_tested_variants) {
  using S = typename Module<T>::Summand_t;
  using B = typename S::Births;
  using D = typename S::Deaths;
  T inf = Module<T>::T_inf;

  const int numParam = 2;
  S s1({1, 6, 2, 5}, {8, 9, 12, 8}, numParam, 0);
  S s2({4, 3}, {inf, inf}, numParam, 2);
  S s3({5, 2, 8, 1}, {13, 7, 14, 6}, numParam, 0);
  S s5(B::inf(numParam), D::inf(numParam), 0);
  S s6(B::minus_inf(numParam), D::minus_inf(numParam), 2);

  Module<T> m;

  m.add_summand(s1);
  m.add_summand(s2);
  m.add_summand(s5);
  m.add_summand(s3);
  m.add_summand(s6);
  m.add_summand(s5);

  char* buffer = new char[1024];
  char* ptr = buffer;

  std::size_t serSize = get_serialization_size_of(m);
  ptr = serialize_value_to_char_buffer(m, ptr);
  BOOST_CHECK_EQUAL(serSize, ptr - buffer);

  Module<T> copy;
  const char* c_ptr = buffer;
  c_ptr = deserialize_value_from_char_buffer(copy, c_ptr);
  BOOST_CHECK_EQUAL(serSize, c_ptr - buffer);
  BOOST_CHECK(m == copy);
}
