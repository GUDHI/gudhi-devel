/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <cstdint>  // for std::uint8_t
#include <cmath>    // std::isnan

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_extended_filtration"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>

struct Low_options : Gudhi::Simplex_tree_options_default {
  typedef float Filtration_value;
  typedef std::uint8_t Vertex_handle;
};

typedef boost::mpl::list<Gudhi::Simplex_tree<>,
                         Gudhi::Simplex_tree<Low_options>> list_of_tested_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(basic_simplex_tree_extend_filtration, Stree, list_of_tested_variants) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "BASIC SIMPLEX TREE EXTEND FILTRATION" << std::endl;
  Stree st;
  
  // Inserted simplex:
  //      5   4
  //      o   o
  //     / \ /
  //    o   o
  //   /2\ /3
  //  o   o
  //  1   0
  st.insert_simplex({0}, 1.0);
  st.insert_simplex({1}, 2.0);
  st.insert_simplex({2}, 3.0);
  st.insert_simplex({3}, 4.0);
  st.insert_simplex({4}, 5.0);
  st.insert_simplex({5}, 6.0);
  st.insert_simplex({0, 2});
  st.insert_simplex({1, 2});
  st.insert_simplex({0, 3});
  st.insert_simplex({2, 5});
  st.insert_simplex({3, 4});
  st.insert_simplex({3, 5});
  Stree copy(st);

  auto efd = st.extend_filtration();

  using Filtration_value = typename Stree::Filtration_value;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({6})), Filtration_value(-3.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({0})), Filtration_value(-2.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({1})), Filtration_value(-1.8));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({2})), Filtration_value(-1.6));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({0, 2})), Filtration_value(-1.6));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({1, 2})), Filtration_value(-1.6));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({3})), Filtration_value(-1.4));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({0, 3})), Filtration_value(-1.4));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({4})), Filtration_value(-1.2));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({3, 4})), Filtration_value(-1.2));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({5})), Filtration_value(-1.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({2, 5})), Filtration_value(-1.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({3, 5})), Filtration_value(-1.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({5, 6})), Filtration_value(1.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({4, 6})), Filtration_value(1.2));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({3, 6})), Filtration_value(1.4));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({3, 4, 6})), Filtration_value(1.4));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({3, 5, 6})), Filtration_value(1.4));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({2, 6})), Filtration_value(1.6));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({2, 5, 6})), Filtration_value(1.6));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({1, 6})), Filtration_value(1.8));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({1, 2, 6})), Filtration_value(1.8));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({0, 6})), Filtration_value(2.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({0, 2, 6})), Filtration_value(2.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({0, 3, 6})), Filtration_value(2.0));

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(efd.minval, Filtration_value(1.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(efd.maxval, Filtration_value(6.0));

  Filtration_value epsilon = 2 * std::numeric_limits<Filtration_value>::epsilon();
  Filtration_value filt;
  Gudhi::Extended_simplex_type est = Gudhi::Extended_simplex_type::EXTRA;
  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({0})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({0})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({1})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({1})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({2})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({2})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({3})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({3})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({4})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({4})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({5})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({5})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({6})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  BOOST_CHECK(std::isnan(filt));
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::EXTRA);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({0, 6})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({0})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::DOWN);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({1, 6})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({1})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::DOWN);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({2, 6})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({2})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::DOWN);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({3, 6})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({3})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::DOWN);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({4, 6})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({4})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::DOWN);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({5, 6})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({5})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::DOWN);

  // Test to the limit - Intervals are [-2., -1.] and [1., 2.]
  std::tie(filt, est) = st.decode_extended_filtration(Filtration_value(-2.1), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  BOOST_CHECK(std::isnan(filt));
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::EXTRA);

  std::tie(filt, est) = st.decode_extended_filtration(Filtration_value(2.1), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  BOOST_CHECK(std::isnan(filt));
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::EXTRA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_extend_same_filtration, Stree, list_of_tested_variants) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "BASIC SIMPLEX TREE EXTEND WITH THE SAME FILTRATION LIMIT CASE" << std::endl;
  Stree st;
  
  // some vertices with the same filtration values
  st.insert_simplex({0}, 1.0);
  st.insert_simplex({1}, 1.0);
  st.insert_simplex({2}, 1.0);
  st.insert_simplex({3}, 1.0);
  st.insert_simplex({4}, 1.0);
  st.insert_simplex({5}, 1.0);
  Stree copy(st);

  auto efd = st.extend_filtration();

  using Filtration_value = typename Stree::Filtration_value;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({0})), Filtration_value(-2.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({1})), Filtration_value(-2.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({2})), Filtration_value(-2.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({3})), Filtration_value(-2.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({4})), Filtration_value(-2.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({5})), Filtration_value(-2.0));
  // Check additional point for coning the simplicial complex
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(st.filtration(st.find({6})), Filtration_value(-3.0));

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(efd.minval, Filtration_value(1.0));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(efd.maxval, Filtration_value(1.0));

  Filtration_value epsilon = 2 * std::numeric_limits<Filtration_value>::epsilon();
  Filtration_value filt;
  Gudhi::Extended_simplex_type est = Gudhi::Extended_simplex_type::EXTRA;
  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({0})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({0})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({1})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({1})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({2})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({2})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({3})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({3})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({4})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({4})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({5})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filt, copy.filtration(copy.find({5})), epsilon);
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::UP);

  std::tie(filt, est) = st.decode_extended_filtration(st.filtration(st.find({6})), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  BOOST_CHECK(std::isnan(filt));
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::EXTRA);

  // Test to the limit - Intervals are [-2., -1.] and [1., 2.]
  std::tie(filt, est) = st.decode_extended_filtration(Filtration_value(-0.9), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  BOOST_CHECK(std::isnan(filt));
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::EXTRA);

  std::tie(filt, est) = st.decode_extended_filtration(Filtration_value(0.9), efd);
  std::clog << filt << " - " << static_cast<int>(est) << std::endl;
  BOOST_CHECK(std::isnan(filt));
  BOOST_CHECK(est == Gudhi::Extended_simplex_type::EXTRA);
}