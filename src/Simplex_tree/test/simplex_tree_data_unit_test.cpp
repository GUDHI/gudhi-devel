/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse, Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - 2024/11 Vincent Rouvreau: Move data test from simplex_tree_unit_test
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <string>
#include <utility>  // for std::move

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_data"
#include <boost/test/unit_test.hpp>

#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

struct Options_with_int_data : Simplex_tree_options_minimal {
  typedef int Simplex_data;
};

BOOST_AUTO_TEST_CASE(simplex_data) {
  Simplex_tree<Options_with_int_data> st;
  st.insert_simplex_and_subfaces({0, 1});
  st.insert_simplex_and_subfaces({2, 1});
  st.insert_simplex_and_subfaces({0, 2});
  st.simplex_data(st.find({0, 1})) = 5;
  st.expansion(3);
  st.simplex_data(st.find({0, 1, 2})) = 4;
  for (auto simplex : st.complex_simplex_range()) {
    std::clog << "(";
    for (auto vertex : st.simplex_vertex_range(simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ")";
    std::clog << " filtration=" << st.filtration(simplex) << " - data='" << st.simplex_data(simplex) << "'\n";
  }
  BOOST_CHECK(st.simplex_data(st.find({0, 1})) == 5);
}

struct Options_with_string_data : Simplex_tree_options_full_featured {
  using Simplex_data = std::string;
};

template<class Simplex_tree_source, class Simplex_tree_copy>
void test_filtration_and_additionnal_data(Simplex_tree_source stree, Simplex_tree_copy stree_copy) {
  for (auto simplex : stree_copy.complex_simplex_range()) {
    std::vector<typename Simplex_tree_source::Vertex_handle> vec_of_vertices {};
    std::clog << "(";
    for (auto vertex : stree_copy.simplex_vertex_range(simplex)) {
      std::clog << vertex << " ";
      vec_of_vertices.push_back(vertex);
    }
    std::clog << ")";
    auto sh = stree.find(vec_of_vertices);
    BOOST_CHECK(sh != stree.null_simplex());
    std::clog << " filtration=" << stree_copy.filtration(simplex) << " versus " << stree.filtration(sh);
    std::clog << " - data='" << stree_copy.simplex_data(simplex) << "' versus '" << stree.simplex_data(sh) << "'\n";
    BOOST_CHECK(stree_copy.filtration(simplex) == stree.filtration(sh));
    BOOST_CHECK(stree_copy.simplex_data(simplex) == stree.simplex_data(sh));
  }
}

BOOST_AUTO_TEST_CASE(simplex_data_copy) {
  Simplex_tree<Options_with_string_data> stree;
  stree.insert_simplex_and_subfaces({0, 1}, 1.);
  stree.simplex_data(stree.find({0, 1})) = "{0, 1}";
  // vertices have a specific algorithm, so let's test one data on a vertex
  stree.simplex_data(stree.find({0})) = "{0}";
  stree.insert_simplex_and_subfaces({2, 1}, 2.);
  stree.simplex_data(stree.find({1, 2})) = "{1, 2}";
  stree.insert_simplex_and_subfaces({0, 2}, 3.);
  stree.simplex_data(stree.find({0, 2})) = "{0, 2}";
  std::clog << "Iterator on stree:\n";
  for (auto simplex : stree.complex_simplex_range()) {
    std::clog << "(";
    for (auto vertex : stree.simplex_vertex_range(simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") filtration=" << stree.filtration(simplex) << " - data='" << stree.simplex_data(simplex) << "'\n";
  }
  {
    Simplex_tree<Options_with_string_data> stree_copy(stree);
    std::clog << "\nIterator on copy constructor:\n";
    test_filtration_and_additionnal_data(stree, stree_copy);
  }
  {
    Simplex_tree<Options_with_string_data> stree_copy = stree;
    std::clog << "\nIterator on copy assignment:\n";
    test_filtration_and_additionnal_data(stree, stree_copy);
  }
  {
    // Needs to copy first for not to lose the original one
    Simplex_tree<Options_with_string_data> stree_copy = stree;
    Simplex_tree<Options_with_string_data> stree_moved(std::move(stree_copy));
    std::clog << "\nIterator on move constructor:\n";
    test_filtration_and_additionnal_data(stree, stree_moved);
  }
  {
    // Needs to copy first for not to lose the original one
    Simplex_tree<Options_with_string_data> stree_copy = stree;
    Simplex_tree<Options_with_string_data> stree_moved = std::move(stree_copy);
    std::clog << "\nIterator on move assignment:\n";
    test_filtration_and_additionnal_data(stree, stree_moved);
  }
  {
    auto trans = [](const typename Simplex_tree<Options_with_string_data>::Filtration_value& fil)
      -> Simplex_tree<Options_with_int_data>::Filtration_value {
        // spoiler: Options_with_int_data does not store filtrations
        if constexpr(Options_with_int_data::store_filtration)
          return fil;
        else
         return {};
      };
    Simplex_tree<Options_with_int_data> stree_copy(stree, trans);
    std::clog << "\nIterator on transformer copy constructor:\n";
    Simplex_tree<Options_with_int_data>::Filtration_value default_filtration_value{};
    Simplex_tree<Options_with_int_data>::Simplex_data default_data_value{};
    for (auto simplex : stree_copy.complex_simplex_range()) {
      std::vector<Simplex_tree<Options_with_int_data>::Vertex_handle> vec_of_vertices {};
      std::clog << "(";
      for (auto vertex : stree_copy.simplex_vertex_range(simplex)) {
        std::clog << vertex << " ";
        vec_of_vertices.push_back(vertex);
      }
      std::clog << ")";
      auto sh = stree.find(vec_of_vertices);
      BOOST_CHECK(sh != stree.null_simplex());
      std::clog << " filtration=" << stree_copy.filtration(simplex) << " versus " << stree.filtration(sh);
      std::clog << " - data='" << stree_copy.simplex_data(simplex) << "' versus '" << stree.simplex_data(sh) << "'\n";
      BOOST_CHECK(stree_copy.filtration(simplex) == default_filtration_value);
      // With transformer constructors, simplex_data is not copied and not initialized, so we cannot test it
      // I just keep it in comment if we find something smarter
      // BOOST_CHECK(stree_copy.simplex_data(simplex) == default_data_value);
    }
  }
}
