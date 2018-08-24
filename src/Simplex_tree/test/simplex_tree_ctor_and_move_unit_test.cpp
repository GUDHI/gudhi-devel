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

  std::cout << "********************************************************************\n";
  std::cout << "* " << msg << "\n";
  std::cout << "* The complex contains " << st.num_simplices() << " simplices";
  std::cout << "   - dimension " << st.dimension() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "   "
              << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) std::cout << "(" << vertex << ")";
    std::cout << std::endl;
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

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF COPY CONSTRUCTOR" << std::endl;

  Simplex_tree st1(st);
  Simplex_tree st2(st);
  print_simplex_filtration(st1, "First copy constructor from the default Simplex_tree");
  print_simplex_filtration(st2, "Second copy constructor from the default Simplex_tree");
  // Cross check
  BOOST_CHECK(st1 == st2);
  BOOST_CHECK(st == st2);
  BOOST_CHECK(st1 == st);

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF COPY ASSIGNMENT" << std::endl;
  Simplex_tree st3;
  st3 = st;
  print_simplex_filtration(st3, "First copy assignment from the default Simplex_tree");
  Simplex_tree st4;
  st4 = st;
  print_simplex_filtration(st4, "Second copy assignment from the default Simplex_tree");

  // Cross check
  BOOST_CHECK(st3 == st4);
  BOOST_CHECK(st == st4);
  BOOST_CHECK(st3 == st);

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF MOVE CONSTRUCTOR" << std::endl;
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

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF MOVE ASSIGNMENT" << std::endl;

  // A swap is a copy ctor of a tmp value, then it uses move assignment
  std::swap(st3, st4);
  print_simplex_filtration(st3, "First move assignment from the default Simplex_tree");
  print_simplex_filtration(st4, "Second move assignment from the default Simplex_tree");
  BOOST_CHECK(st3 == st4);
  BOOST_CHECK(st == st4);
  BOOST_CHECK(st3 == st);

}
