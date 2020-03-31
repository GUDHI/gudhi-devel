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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_constructor_and_move"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

typedef boost::mpl::list<Simplex_tree<>, Simplex_tree<Simplex_tree_options_fast_persistence>> list_of_tested_variants;

template<typename Simplex_tree>
void print_simplex_filtration(Simplex_tree& st, const std::string& msg) {
  // Required before browsing through filtration values
  st.initialize_filtration();

  std::clog << "********************************************************************\n";
  std::clog << "* " << msg << "\n";
  std::clog << "* The complex contains " << st.num_simplices() << " simplices";
  std::clog << "   - dimension " << st.dimension() << "\n";
  std::clog << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::clog << "   "
              << "[" << st.filtration(f_simplex) << "] ";
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

  st = st;
  print_simplex_filtration(st4, "Third self copy assignment from the default Simplex_tree");

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

  st = std::move(st);
  print_simplex_filtration(st, "Third self move assignment from the default Simplex_tree");

  BOOST_CHECK(st7 == st);

}
