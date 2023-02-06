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
#include <cstring>  // for std::size_t and strncmp
#include <random>
#include <iterator>  // for std::distance
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_serialization"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Simplex_tree/serialization_utils.h>  // for Gudhi::simplex_tree::get_serialization_size, serialize, deserialize
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

  std::vector<char> buffer;
  // 3 simplices ({0}, {1}, {2}) and its filtration values
  serialize(static_cast<Vertex_type>(3)   , buffer);
  serialize(static_cast<Vertex_type>(0)   , buffer);
  serialize(st.filtration(st.find({0}))   , buffer);
  serialize(static_cast<Vertex_type>(1)   , buffer);
  serialize(st.filtration(st.find({1}))   , buffer);
  serialize(static_cast<Vertex_type>(2)   , buffer);
  serialize(st.filtration(st.find({2}))   , buffer);
  // 1 simplex (2) from {0, 2} and its filtration values
  serialize(static_cast<Vertex_type>(1)   , buffer);
  serialize(static_cast<Vertex_type>(2)   , buffer);
  serialize(st.filtration(st.find({0, 2})), buffer);
  serialize(static_cast<Vertex_type>(0)   , buffer);  // (0, 2) end of leaf
  serialize(static_cast<Vertex_type>(0)   , buffer);  // (1) end of leaf
  serialize(static_cast<Vertex_type>(0)   , buffer);  // (2) end of leaf

  std::clog << "Serialization size in bytes = " << buffer.size() << std::endl;
  // Sizes are expressed in bytes
  const std::size_t vertex_size = sizeof(Vertex_type);
  const std::size_t filtration_size = sizeof(Filtration_type);
  const std::size_t serialization_size = vertex_size + st.num_simplices() * (2 * vertex_size + filtration_size);
  BOOST_CHECK(serialization_size == buffer.size());

  Vertex_type vertex = 0;
  Filtration_type filtration = 0;
  // Reset position pointer at start
  auto position_ptr = std::cbegin(buffer);
  // 3 simplices ({0}, {1}, {2}) and its filtration values
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 3);
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 0);
  position_ptr = deserialize(position_ptr, filtration);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({0})));
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 1);
  position_ptr = deserialize(position_ptr, filtration);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({1})));
  position_ptr = deserialize(position_ptr, vertex);
  BOOST_CHECK(vertex == 2);
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

  std::size_t index = std::distance(std::cbegin(buffer), position_ptr);
  std::clog << "Deserialization size in bytes = " << serialization_size << " - index = " << index << std::endl;
  BOOST_CHECK(serialization_size == index);
  BOOST_CHECK(position_ptr == std::cend(buffer));

  std::vector<char> stree_buffer = st.serialize();

  std::clog << "Serialization (from the simplex tree) size in bytes = " << stree_buffer.size() << std::endl;
  BOOST_CHECK(serialization_size == stree_buffer.size());

  std::clog << "\nSerialization (from the simplex tree):\n";
  for (std::size_t idx = 0; idx < serialization_size; idx++) {
    std::clog << std::setfill('0') 
              << std::setw(2) 
              << std::uppercase 
              << std::hex << (0xFF & stree_buffer[idx]) << " " << std::dec;
  }
  std::clog << "\nSerialization (manual):\n";
  for (std::size_t idx = 0; idx < serialization_size; idx++) {
    std::clog << std::setfill('0') 
              << std::setw(2) 
              << std::uppercase 
              << std::hex << (0xFF & buffer[idx]) << " " << std::dec;
  }
  BOOST_CHECK(stree_buffer == buffer);

  Stree* st_from_buffer = Stree::deserialize(stree_buffer);
  std::clog << std::endl << std::endl << "Iterator on simplices:\n";
  for (auto simplex : st_from_buffer->complex_simplex_range()) {
    std::clog << "   ";
    for (auto vertex : st_from_buffer->simplex_vertex_range(simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << " - filtration = " << st_from_buffer->filtration(simplex) << std::endl;
  }
  BOOST_CHECK(*st_from_buffer == st);
}
