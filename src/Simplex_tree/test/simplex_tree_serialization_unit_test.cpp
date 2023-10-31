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
#include <cstdint>  // for std::uint8_t
#include <iomanip>  // for std::setfill, setw
#include <ios>  // for std::hex, uppercase

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_serialization"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Simplex_tree/serialization_utils.h>  // for de/serialize_trivial
#include <gudhi/Unitary_tests_utils.h>  // for GUDHI_TEST_FLOAT_EQUALITY_CHECK

using namespace Gudhi;
using namespace Gudhi::simplex_tree;

struct Low_options : Gudhi::Simplex_tree_options_full_featured {
  // Implicitly use 0 as filtration value for all simplices
  static const bool store_filtration = false;
  // The persistence algorithm needs this
  static const bool store_key = true;
  // I have few vertices
  typedef std::uint8_t Vertex_handle;
  // Maximum number of simplices to compute persistence is 2^8 - 1 = 255. One is reserved for null_key
  typedef std::uint8_t Simplex_key;
};

struct Stable_options : Gudhi::Simplex_tree_options_full_featured {
  //disabled by default.
  static const bool stable_simplex_handles = true;
};

typedef boost::mpl::list<Simplex_tree<>,
                         Simplex_tree<Simplex_tree_options_fast_persistence>,
                         Simplex_tree<Low_options>,
                         Simplex_tree<Simplex_tree_options_fast_cofaces>,
                         Simplex_tree<Stable_options> > list_of_tested_variants;

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

  if (Stree::Options::store_filtration)
    st.insert_simplex({0},    random_filtration<Filtration_type>());
  else
    st.insert_simplex({0});
  if (Stree::Options::store_filtration)
    st.insert_simplex({1},    random_filtration<Filtration_type>());
  else
    st.insert_simplex({1});
  if (Stree::Options::store_filtration)
    st.insert_simplex({2},    random_filtration<Filtration_type>());
  else
    st.insert_simplex({2});
  if (Stree::Options::store_filtration)
    st.insert_simplex({0, 2}, random_filtration<Filtration_type>());
  else
    st.insert_simplex({0, 2});

  char* buffer = new char[256];
  char* ptr = buffer;
  // 3 simplices ({0}, {1}, {2}) and its filtration values
  ptr = serialize_trivial(static_cast<Vertex_type>(3)     , ptr);
  ptr = serialize_trivial(static_cast<Vertex_type>(0)     , ptr);
  if (Stree::Options::store_filtration)
    ptr = serialize_trivial(st.filtration(st.find({0}))   , ptr);
  ptr = serialize_trivial(static_cast<Vertex_type>(1)     , ptr);
  if (Stree::Options::store_filtration)
    ptr = serialize_trivial(st.filtration(st.find({1}))   , ptr);
  ptr = serialize_trivial(static_cast<Vertex_type>(2)     , ptr);
  if (Stree::Options::store_filtration)
    ptr = serialize_trivial(st.filtration(st.find({2}))   , ptr);
  // 1 simplex (2) from {0, 2} and its filtration values
  ptr = serialize_trivial(static_cast<Vertex_type>(1)     , ptr);
  ptr = serialize_trivial(static_cast<Vertex_type>(2)     , ptr);
  if (Stree::Options::store_filtration)
    ptr = serialize_trivial(st.filtration(st.find({0, 2})), ptr);
  ptr = serialize_trivial(static_cast<Vertex_type>(0)     , ptr);  // (0, 2) end of leaf
  ptr = serialize_trivial(static_cast<Vertex_type>(0)     , ptr);  // (1) end of leaf
  ptr = serialize_trivial(static_cast<Vertex_type>(0)     , ptr);  // (2) end of leaf

  const std::size_t buffer_size = (ptr - buffer);
  std::clog << "Serialization size in bytes = " << buffer_size << std::endl;
  // Sizes are expressed in bytes
  const std::size_t vertex_size = sizeof(Vertex_type);
  const std::size_t filtration_size = Stree::Options::store_filtration ? sizeof(Filtration_type) : 0;
  const std::size_t serialization_size = vertex_size + st.num_simplices() * (2 * vertex_size + filtration_size);
  BOOST_CHECK(serialization_size == buffer_size);

  Vertex_type vertex = 0;
  Filtration_type filtration = 0;
  // Reset position pointer at start
  const char* c_ptr = buffer;
  // 3 simplices ({0}, {1}, {2}) and its filtration values
  c_ptr = deserialize_trivial(vertex, c_ptr);
  BOOST_CHECK(vertex == 3);
  c_ptr = deserialize_trivial(vertex, c_ptr);
  BOOST_CHECK(vertex == 0);
  if (Stree::Options::store_filtration) {
    c_ptr = deserialize_trivial(filtration, c_ptr);
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({0})));
  }
  c_ptr = deserialize_trivial(vertex, c_ptr);
  BOOST_CHECK(vertex == 1);
  if (Stree::Options::store_filtration) {
    c_ptr = deserialize_trivial(filtration, c_ptr);
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({1})));
  }
  c_ptr = deserialize_trivial(vertex, c_ptr);
  BOOST_CHECK(vertex == 2);
  if (Stree::Options::store_filtration) {
    c_ptr = deserialize_trivial(filtration, c_ptr);
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({2})));
  }
  // 1 simplex (2) from {0, 2} and its filtration values
  c_ptr = deserialize_trivial(vertex, c_ptr);
  BOOST_CHECK(vertex == 1);
  c_ptr = deserialize_trivial(vertex, c_ptr);
  BOOST_CHECK(vertex == 2);
  if (Stree::Options::store_filtration) {
    c_ptr = deserialize_trivial(filtration, c_ptr);
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(filtration, st.filtration(st.find({0, 2})));
  }
  c_ptr = deserialize_trivial(vertex, c_ptr);  // (0, 2) end of leaf
  BOOST_CHECK(vertex == 0);
  c_ptr = deserialize_trivial(vertex, c_ptr);  // (1) end of leaf
  BOOST_CHECK(vertex == 0);
  c_ptr = deserialize_trivial(vertex, c_ptr);  // (2) end of leaf
  BOOST_CHECK(vertex == 0);

  std::size_t index = static_cast<std::size_t>(c_ptr - buffer);
  std::clog << "Deserialization size in bytes = " << serialization_size << " - index = " << index << std::endl;
  BOOST_CHECK(serialization_size == index);

  const std::size_t stree_buffer_size = st.get_serialization_size();
  std::clog << "Serialization (from the simplex tree) size in bytes = " << stree_buffer_size << std::endl;
  BOOST_CHECK(serialization_size == stree_buffer_size);
  char* stree_buffer = new char[stree_buffer_size];
  st.serialize(stree_buffer, stree_buffer_size);


  std::clog << "\nSerialization (from the simplex tree):\n";
  for (std::size_t idx = 0; idx < stree_buffer_size; idx++) {
    std::clog << std::setfill('0') 
              << std::setw(2) 
              << std::uppercase 
              << std::hex << (0xFF & stree_buffer[idx]) << " " << std::dec;
  }
  std::clog << "\nSerialization (manual):\n";
  for (std::size_t idx = 0; idx < stree_buffer_size; idx++) {
    std::clog << std::setfill('0') 
              << std::setw(2) 
              << std::uppercase 
              << std::hex << (0xFF & buffer[idx]) << " " << std::dec;
  }
  BOOST_CHECK(strncmp(stree_buffer, buffer, stree_buffer_size) == 0);

  Stree st_from_buffer;
  st_from_buffer.deserialize(stree_buffer, serialization_size);
  std::clog << std::endl << std::endl << "Iterator on simplices:\n";
  for (auto simplex : st_from_buffer.complex_simplex_range()) {
    std::clog << "   ";
    for (auto vertex : st_from_buffer.simplex_vertex_range(simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << " - filtration = " << st_from_buffer.filtration(simplex) << std::endl;
  }
  BOOST_CHECK(st_from_buffer == st);

  delete[] buffer;
  delete[] stree_buffer;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_serialization_exception, Stree, list_of_tested_variants) {
  const std::size_t too_long_buffer_size = 256;
  char* buffer = new char[too_long_buffer_size];

  std::clog << "Serialization of a too long buffer:\n";
  for (std::size_t idx = 0; idx < too_long_buffer_size; idx++) {
    std::clog << std::setfill('0') 
              << std::setw(2) 
              << std::uppercase 
              << std::hex << (0xFF & buffer[idx]) << " " << std::dec;
  }
  std::clog << std::endl;
  Stree st;
  BOOST_CHECK_THROW(st.serialize(buffer, too_long_buffer_size), std::invalid_argument);
  delete[] buffer;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_deserialization_exception, Stree, list_of_tested_variants) {
  const std::size_t too_long_buffer_size = 256;
  char buffer[too_long_buffer_size]{};

  std::clog << "Deserialization of a too long buffer:\n";
  for (std::size_t idx = 0; idx < too_long_buffer_size; idx++) {
    std::clog << std::setfill('0') 
              << std::setw(2) 
              << std::uppercase 
              << std::hex << (0xFF & buffer[idx]) << " " << std::dec;
  }
  std::clog << std::endl;
  Stree st;
  BOOST_CHECK_THROW(st.deserialize(buffer, too_long_buffer_size), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_empty_serialize_deserialize, Stree, list_of_tested_variants) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "EMPTY SIMPLEX TREE SERIALIZATION/DESERIALIZATION" << std::endl;
  Stree st;
  std::clog << "Empty Simplex_tree dimension = " << st.dimension() << std::endl;

  const std::size_t stree_buffer_size = st.get_serialization_size();
  std::clog << "Serialization (from the simplex tree) size in bytes = " << stree_buffer_size << std::endl;
  char* stree_buffer = new char[stree_buffer_size];

  st.serialize(stree_buffer, stree_buffer_size);

  Stree st_from_buffer;
  st_from_buffer.deserialize(stree_buffer, stree_buffer_size);
  std::clog << "Empty Simplex_tree dimension = " << st_from_buffer.dimension() << std::endl;

  BOOST_CHECK(st_from_buffer == st);

  delete[] stree_buffer;
}


#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_non_empty_deserialize_throw, Stree, list_of_tested_variants) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "NON EMPTY SIMPLEX TREE DESERIALIZATION EXCEPTION" << std::endl;
  Stree st;

  const std::size_t stree_buffer_size = st.get_serialization_size();
  std::clog << "Serialization (from the simplex tree) size in bytes = " << stree_buffer_size << std::endl;
  char* stree_buffer = new char[stree_buffer_size];

  st.serialize(stree_buffer, stree_buffer_size);

  Stree st_from_buffer;
  if (Stree::Options::store_filtration)
    st_from_buffer.insert_simplex({0}, random_filtration<typename Stree::Filtration_value>());
  else
    st_from_buffer.insert_simplex({0});
  std::clog << "Check exception throw in debug mode" << std::endl;
  // throw exception because st_from_buffer is not empty
  BOOST_CHECK_THROW (st_from_buffer.deserialize(stree_buffer, stree_buffer_size),
                     std::logic_error);
  delete[] stree_buffer;
}
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_serialization_and_cofaces, Stree, list_of_tested_variants) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "BASIC SIMPLEX TREE SERIALIZATION/DESERIALIZATION" << std::endl;
  Stree st;

  st.insert_simplex_and_subfaces({2, 1, 0});
  st.insert_simplex_and_subfaces({3, 0});
  st.insert_simplex_and_subfaces({3, 4, 5});
  st.insert_simplex_and_subfaces({0, 1, 6, 7});
  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o           */
  /*  2   0   3\X/4           */
  /*            o             */
  /*            5             */

  const std::size_t stree_buffer_size = st.get_serialization_size();
  std::clog << "Serialization (from the simplex tree) size in bytes = " << stree_buffer_size << std::endl;
  char* stree_buffer = new char[stree_buffer_size];

  st.serialize(stree_buffer, stree_buffer_size);

  Stree st_from_buffer;
  st_from_buffer.deserialize(stree_buffer, stree_buffer_size);
  delete[] stree_buffer;

  int num_stars = 0;
  for (auto coface : st_from_buffer.star_simplex_range(st.find({1, 0}))) {
    std::clog << "coface";
    for (auto vertex : st_from_buffer.simplex_vertex_range(coface)) {
      std::clog << " " << vertex;
    }
    std::clog << "\n";
    num_stars++;
  }
  // [([0, 1], 0.0), ([0, 1, 2], 0.0), ([0, 1, 6], 0.0), ([0, 1, 6, 7], 0.0), ([0, 1, 7], 0.0)]
  BOOST_CHECK(num_stars == 5);

}
