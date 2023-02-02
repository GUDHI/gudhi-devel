/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <cstddef>  // for std::size_t
#include <random>
#include <iterator>  // for std::distance

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_serialization"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "gudhi/Simplex_tree.h"
#include "gudhi/Simplex_tree/serialization_utils.h"  // for Gudhi::simplex_tree::get_serialization_size, serialize, deserialize
#include <gudhi/Unitary_tests_utils.h>  // for GUDHI_TEST_FLOAT_EQUALITY_CHECK

using namespace Gudhi;
using namespace Gudhi::simplex_tree;

typedef boost::mpl::list<Simplex_tree<>, Simplex_tree<Simplex_tree_options_fast_persistence>> list_of_tested_variants;

template<class Filtration_type>
Filtration_type random_filtration(Filtration_type lower_bound = 0, Filtration_type upper_bound = 1) {
  std::uniform_real_distribution<Filtration_type> unif(lower_bound, upper_bound);
  std::random_device rand_dev;
  std::mt19937 rand_engine(rand_dev());
  
  return unif(rand_engine);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(basic_simplex_tree_serialization, Stree, list_of_tested_variants) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "BASIC SIMPLEX TREE SERIALIZATION/DESERIALIZATION" << std::endl;
  Stree st;

  using Filtration_type = typename Stree::Filtration_value;
  using Vertex_type = typename Stree::Vertex_handle;

  st.insert_simplex({0},    random_filtration<Filtration_type>());
  st.insert_simplex({1},    random_filtration<Filtration_type>());
  st.insert_simplex({2},    random_filtration<Filtration_type>());
  st.insert_simplex({0, 2}, random_filtration<Filtration_type>());

  std::clog << "Number of simplices = " << st.num_simplices() << std::endl;
  const std::size_t serialization_size_in_bytes = get_serialization_size<Stree>(st.num_simplices());

  char* serial = new char[serialization_size_in_bytes];
  // Set position pointer at start
  char* position_ptr = serial;
  // 3 simplices ({0}, {1}, {2}) and its filtration values
  position_ptr = serialize(position_ptr, static_cast<Vertex_type>(3));
  position_ptr = serialize(position_ptr, static_cast<Vertex_type>(0));
  position_ptr = serialize(position_ptr, static_cast<Vertex_type>(1));
  position_ptr = serialize(position_ptr, static_cast<Vertex_type>(2));
  position_ptr = serialize(position_ptr, st.filtration(st.find({0})));
  position_ptr = serialize(position_ptr, st.filtration(st.find({1})));
  position_ptr = serialize(position_ptr, st.filtration(st.find({2})));
  // 1 simplex (2) from {0, 2} and its filtration values
  position_ptr = serialize(position_ptr, static_cast<Vertex_type>(1));
  position_ptr = serialize(position_ptr, static_cast<Vertex_type>(2));
  position_ptr = serialize(position_ptr, st.filtration(st.find({0, 2})));
  position_ptr = serialize(position_ptr, static_cast<Vertex_type>(0));  // (0, 2) end of leaf
  position_ptr = serialize(position_ptr, static_cast<Vertex_type>(0));  // (1) end of leaf
  position_ptr = serialize(position_ptr, static_cast<Vertex_type>(0));  // (2) end of leaf

  std::size_t index = std::distance(serial, position_ptr);
  std::clog << "Serialization size in bytes = " << serialization_size_in_bytes << " - index = " << index << std::endl;
  BOOST_CHECK(serialization_size_in_bytes == index);

  Vertex_type vertex = 0;
  Filtration_type filtration = 0;
  // Reset position pointer at start
  position_ptr = serial;
  // 3 simplices ({0}, {1}, {2}) and its filtration values
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 3);
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 0);
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 1);
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 2);
  position_ptr = deserialize(position_ptr, filtration);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({0})));
  position_ptr = deserialize(position_ptr, filtration);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({1})));
  position_ptr = deserialize(position_ptr, filtration);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({2})));
  // 1 simplex (2) from {0, 2} and its filtration values
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 1);
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 2);
  position_ptr = deserialize(position_ptr, filtration);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({0, 2})));
  position_ptr = deserialize(position_ptr, vertex);  // (0, 2) end of leaf
  BOOST_CHECK(vertex == 0);
  position_ptr = deserialize(position_ptr, vertex);  // (1) end of leaf
  BOOST_CHECK(vertex == 0);
  position_ptr = deserialize(position_ptr, vertex);  // (2) end of leaf
  BOOST_CHECK(vertex == 0);

  index = std::distance(serial, position_ptr);
  std::clog << "Deserialization size in bytes = " << serialization_size_in_bytes << " - index = " << index << std::endl;
  BOOST_CHECK(serialization_size_in_bytes == index);

  delete[] serial;
}
