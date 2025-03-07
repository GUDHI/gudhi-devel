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
#include <vector>
#include <string>
#include <iterator>  // for std::distance
#include <algorithm>  // for std::equal
#include <utility>  // for std::move

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_constructor_and_move"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/Simplex_tree.h"

#include "test_vector_filtration_simplex_tree.h"

using namespace Gudhi;

typedef boost::mpl::list<Simplex_tree<>,
                         Simplex_tree<Simplex_tree_options_fast_persistence>,
                         Simplex_tree<Simplex_tree_options_full_featured> > list_of_tested_variants;

std::string print_filtration_value(double fil){
  return std::to_string(fil);
}

std::string print_filtration_value(std::vector<int> fil){
  std::string ss;
  for (auto val : fil) ss += std::to_string(val) + " ";
  ss.pop_back();
  return ss;
}

template<typename Simplex_tree>
void print_simplex_filtration(const Simplex_tree& st, const std::string& msg) {
  // Required before browsing through filtration values
  st.initialize_filtration();

  std::clog << "********************************************************************\n";
  std::clog << "* " << msg << "\n";
  std::clog << "* The complex contains " << st.num_simplices() << " simplices";
  std::clog << "   - dimension " << st.dimension() << "\n";
  std::clog << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::clog << "   "
              << "[" << print_filtration_value(st.filtration(f_simplex)) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_copy_constructor, Simplex_tree, list_of_tested_variants) {
  Simplex_tree st;

  st.insert_simplex_and_subfaces({2, 1, 0}, 3.0);
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, 4.0);
  st.insert_simplex_and_subfaces({3, 0}, 2.0);
  st.insert_simplex_and_subfaces({3, 4, 5}, 3.0);
  st.insert_simplex_and_subfaces({8}, 1.0);
  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o   o       */
  /*  2   0   3\X/4   8       */
  /*            o             */
  /*            5             */
  /*                          */
  /* In other words:          */
  /*   A facet  [2,1,0]       */
  /*   An edge  [0,3]         */
  /*   A facet  [3,4,5]       */
  /*   A cell   [0,1,6,7]     */
  /*   A vertex [8]           */

  print_simplex_filtration(st, "Default Simplex_tree is initialized");

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF COPY CONSTRUCTOR" << std::endl;

  Simplex_tree st1(st);
  Simplex_tree st2(st);
  print_simplex_filtration(st1, "First copy constructor from the default Simplex_tree");
  print_simplex_filtration(st2, "Second copy constructor from the default Simplex_tree");
  // Cross check
  BOOST_CHECK(st1 == st2);
  BOOST_CHECK(st == st2);
  BOOST_CHECK(st1 == st);

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF COPY ASSIGNMENT" << std::endl;
  Simplex_tree st3;
  // To check there is no memory leak
  st3.insert_simplex_and_subfaces({9, 10, 11}, 200.0);
  st3 = st;
  print_simplex_filtration(st3, "First copy assignment from the default Simplex_tree");
  Simplex_tree st4;
  st4 = st;
  print_simplex_filtration(st4, "Second copy assignment from the default Simplex_tree");

  // Cross check
  BOOST_CHECK(st3 == st4);
  BOOST_CHECK(st == st4);
  BOOST_CHECK(st3 == st);

#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
  st = st;
#ifdef __clang__
#pragma GCC diagnostic pop
#endif

  print_simplex_filtration(st, "Third self copy assignment from the default Simplex_tree");

  BOOST_CHECK(st3 == st);

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF MOVE CONSTRUCTOR" << std::endl;
  Simplex_tree st5(std::move(st1));
  print_simplex_filtration(st5, "First move constructor from the default Simplex_tree");
  print_simplex_filtration(st1, "First moved Simplex_tree shall be empty");
  Simplex_tree st6(std::move(st2));
  print_simplex_filtration(st6, "Second move constructor from the default Simplex_tree");
  print_simplex_filtration(st2, "Second moved Simplex_tree shall be empty");

  // Cross check
  BOOST_CHECK(st5 == st6);
  BOOST_CHECK(st == st6);
  BOOST_CHECK(st5 == st);

  Simplex_tree empty_st;
  BOOST_CHECK(st1 == st2);
  BOOST_CHECK(empty_st == st2);
  BOOST_CHECK(st1 == empty_st);

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF MOVE ASSIGNMENT" << std::endl;

  Simplex_tree st7;
  // To check there is no memory leak
  st7.insert_simplex_and_subfaces({9, 10, 11}, 200.0);
  st7 = std::move(st3);
  print_simplex_filtration(st7, "First move assignment from the default Simplex_tree");
  Simplex_tree st8;
  st8 = std::move(st4);
  print_simplex_filtration(st8, "Second move assignment from the default Simplex_tree");

  // Cross check
  BOOST_CHECK(st7 == st8);
  BOOST_CHECK(st == st8);
  BOOST_CHECK(st7 == st);

#if __GNUC__ >= 13
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-move"
#endif
  st = std::move(st);
#if __GNUC__ >= 13
#pragma GCC diagnostic pop
#endif
  print_simplex_filtration(st, "Third self move assignment from the default Simplex_tree");

  BOOST_CHECK(st7 == st);

}

typedef boost::mpl::list<Simplex_tree<Simplex_tree_options_custom_fil_values_default>,
                         Simplex_tree<Simplex_tree_options_custom_fil_values_fast_persistence>,
                         Simplex_tree<Simplex_tree_options_custom_fil_values_full_featured> >
    list_of_custom_fil_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_custom_copy_constructor, typeST, list_of_custom_fil_variants) {
  typeST st;
  Simplex_tree<> st_witness;

  st.insert_simplex_and_subfaces({2, 1, 0}, {3, 1});
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, {4, 1});
  st.insert_simplex_and_subfaces({3, 0}, {2, 1});
  st.insert_simplex_and_subfaces({3, 4, 5}, {3, 1});
  st.insert_simplex_and_subfaces({8}, {1, 1});

  st_witness.insert_simplex_and_subfaces({2, 1, 0}, 3.0);
  st_witness.insert_simplex_and_subfaces({0, 1, 6, 7}, 4.0);
  st_witness.insert_simplex_and_subfaces({3, 0}, 2.0);
  st_witness.insert_simplex_and_subfaces({3, 4, 5}, 3.0);
  st_witness.insert_simplex_and_subfaces({8}, 1.0);
  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o   o       */
  /*  2   0   3\X/4   8       */
  /*            o             */
  /*            5             */
  /*                          */
  /* In other words:          */
  /*   A facet  [2,1,0]       */
  /*   An edge  [0,3]         */
  /*   A facet  [3,4,5]       */
  /*   A cell   [0,1,6,7]     */
  /*   A vertex [8]           */

  print_simplex_filtration(st, "Simplex_tree is initialized");

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF CUSTOM COPY CONSTRUCTOR" << std::endl;

  auto trans = [](const typename typeST::Filtration_value& fil) -> Simplex_tree<>::Filtration_value {
    return fil[0];
  };

  Simplex_tree<> st1(st, trans);
  Simplex_tree<> st2(st, trans);
  print_simplex_filtration(st1, "First custom copy constructor from the Simplex_tree");
  print_simplex_filtration(st2, "Second custom copy constructor from the Simplex_tree");
  // Cross check
  BOOST_CHECK(st1 == st2);
  BOOST_CHECK(st_witness == st2);
  BOOST_CHECK(st1 == st_witness);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_custom_copy_constructor_key, typeST, list_of_custom_fil_variants)
{
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF CUSTOM COPY CONSTRUCTOR WITH KEY VALUES" << std::endl;

  typeST st;

  st.insert_simplex_and_subfaces({2, 1, 0}, {3, 1});
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, {4, 1});
  st.insert_simplex_and_subfaces({3, 0}, {2, 1});
  st.insert_simplex_and_subfaces({3, 4, 5}, {3, 1});
  st.insert_simplex_and_subfaces({8}, {1, 1});

  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o   o       */
  /*  2   0   3\X/4   8       */
  /*            o             */
  /*            5             */
  /*                          */
  /* In other words:          */
  /*   A facet  [2,1,0]       */
  /*   An edge  [0,3]         */
  /*   A facet  [3,4,5]       */
  /*   A cell   [0,1,6,7]     */
  /*   A vertex [8]           */

  if constexpr (typeST::Options::store_key){
    for (auto f_simplex : st.complex_simplex_range()) {
      std::int32_t key = 1;
      for (auto filt : st.filtration(f_simplex)) {
        key *= filt;
      }
      st.assign_key(f_simplex, key);
    }
  }

  auto trans = [](const typename typeST::Filtration_value& fil)
      -> Simplex_tree<Simplex_tree_options_custom_fil_values_full_featured>::Filtration_value {
    Simplex_tree<Simplex_tree_options_custom_fil_values_full_featured>::Filtration_value copy(std::begin(fil),
                                                                                              std::end(fil));
    std::sort(copy.begin(), copy.end());
    return copy;
  };

  Simplex_tree<Simplex_tree_options_custom_fil_values_full_featured> st_trans(st, trans);
  for (auto f_simplex : st_trans.complex_simplex_range()) {
    auto filtrations = st_trans.filtration(f_simplex);
    BOOST_CHECK(std::is_sorted(std::begin(filtrations), std::end(filtrations)));

    if constexpr (typeST::Options::store_key){
      std::int32_t key = 1;
      for (auto filt : st_trans.filtration(f_simplex)) {
        key *= filt;
      }
      std::clog << "key = " << st_trans.key(f_simplex) << " vs " << key << std::endl;

      BOOST_CHECK(st_trans.key(f_simplex) == key);
    } else {
      std::clog << "No key stored initially: key = " << st_trans.key(f_simplex) << std::endl;

      BOOST_CHECK(st_trans.key(f_simplex) == -1);
    }
  }
}

template<typename Simplex_tree>
std::vector<std::vector<typename Simplex_tree::Vertex_handle>> get_star(const Simplex_tree& st) {
  std::vector<std::vector<typename Simplex_tree::Vertex_handle>> output;
  auto sh = st.find({0,1});
  if (sh != st.null_simplex())
  {
    auto stars = st.star_simplex_range(sh);
    output.resize(std::distance(stars.begin(), stars.end()));
    for (auto simplex = stars.begin(); simplex != stars.end(); ++simplex) {
      typename Simplex_tree::Simplex_vertex_range rg = st.simplex_vertex_range(*simplex);
      std::clog << "(";
      for (auto vertex = rg.begin(); vertex != rg.end(); ++vertex) {
        std::clog << *vertex << ", ";
      }
      std::clog << ")" << std::endl;
      output.emplace_back(rg.begin(), rg.end());
    }
  }
  return output;
}

BOOST_AUTO_TEST_CASE(simplex_fast_cofaces_rule_of_five) {
  // Only for fast cofaces version to check the data structure for this feature is up to date
  using STree = Simplex_tree<Simplex_tree_options_full_featured>;
  STree st;

  st.insert_simplex_and_subfaces({2, 1, 0}, 3.0);
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, 4.0);
  st.insert_simplex_and_subfaces({3, 0}, 2.0);
  st.insert_simplex_and_subfaces({3, 4, 5}, 3.0);
  st.insert_simplex_and_subfaces({8}, 1.0);
  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o   o       */
  /*  2   0   3\X/4   8       */
  /*            o             */
  /*            5             */
  /*                          */
  /* In other words:          */
  /*   A facet  [2,1,0]       */
  /*   An edge  [0,3]         */
  /*   A facet  [3,4,5]       */
  /*   A cell   [0,1,6,7]     */
  /*   A vertex [8]           */

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF COPY CONSTRUCTOR FOR COFACES" << std::endl;

  STree st1(st);

  std::clog << "get_star on the original simplex" << std::endl;
  auto stars = get_star(st);
  std::clog << "get_star on the copy by construction" << std::endl;
  auto stars1 = get_star(st1);
  BOOST_CHECK(stars.size() == stars1.size());
  BOOST_CHECK(std::equal(stars.begin(), stars.begin() + stars.size(), stars1.begin()));

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF COPY ASSIGNMENT FOR COFACES" << std::endl;
  STree st2;
  // To check there is no memory leak
  st2.insert_simplex_and_subfaces({9, 10, 11}, 200.0);
  st2 = st;

  std::clog << "get_star on the copy by assignment" << std::endl;
  auto stars2 = get_star(st2);

  BOOST_CHECK(stars.size() == stars2.size());
  BOOST_CHECK(std::equal(stars.begin(), stars.begin() + stars.size(), stars2.begin()));

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF MOVE CONSTRUCTOR FOR COFACES" << std::endl;
  STree st3(std::move(st1));

  std::clog << "get_star on the move by construction" << std::endl;
  auto stars3 = get_star(st3);
  std::clog << "what about the moved one ?" << std::endl;
  auto stars1_moved = get_star(st1);

  BOOST_CHECK(stars.size() == stars3.size());
  BOOST_CHECK(std::equal(stars.begin(), stars.begin() + stars.size(), stars3.begin()));
  BOOST_CHECK(stars1_moved.size() == 0);

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF MOVE ASSIGNMENT FOR COFACES" << std::endl;

  STree st4;
  // To check there is no memory leak
  st4.insert_simplex_and_subfaces({9, 10, 11}, 200.0);
  st4 = std::move(st2);

  std::clog << "get_star on the move by assignment" << std::endl;
  auto stars4 = get_star(st4);
  std::clog << "what about the moved one ?" << std::endl;
  auto stars2_moved = get_star(st2);

  BOOST_CHECK(stars.size() == stars4.size());
  BOOST_CHECK(std::equal(stars.begin(), stars.begin() + stars.size(), stars4.begin()));
  BOOST_CHECK(stars2_moved.size() == 0);
}
