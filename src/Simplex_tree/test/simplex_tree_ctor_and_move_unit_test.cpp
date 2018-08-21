#include <iostream>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_constructor_and_move"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

typedef boost::mpl::list<Simplex_tree<>, Simplex_tree<Simplex_tree_options_fast_persistence>> list_of_tested_variants;

template<class SimplicialComplex>
SimplicialComplex move_it(SimplicialComplex sc) {
  return sc;
}


BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_copy_constructor, Simplex_tree, list_of_tested_variants) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF COPY CONSTRUCTOR" << std::endl;
  Simplex_tree st;

  st.insert_simplex_and_subfaces({2, 1, 0}, 3.0);
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, 4.0);
  st.insert_simplex_and_subfaces({3, 0}, 2.0);
  st.insert_simplex_and_subfaces({3, 4, 5}, 3.0);
  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o           */
  /*  2   0   3\X/4           */
  /*            o             */
  /*            5             */
  /*                          */
  /* In other words:          */
  /*   A facet [2,1,0]        */
  /*   An edge [0,3]          */
  /*   A facet [3,4,5]        */
  /*   A cell  [0,1,6,7]      */

  Simplex_tree st1(st);
  Simplex_tree st2(st);
  // Cross check
  BOOST_CHECK(st1 == st2);
  BOOST_CHECK(st == st2);
  BOOST_CHECK(st1 == st);

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF COPY ASSIGNMENT" << std::endl;
  Simplex_tree st3;
  st3 = st;
  Simplex_tree st4;
  st4 = st;
  // Cross check
  BOOST_CHECK(st3 == st4);
  BOOST_CHECK(st == st4);
  BOOST_CHECK(st3 == st);

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF MOVE CONSTRUCTOR" << std::endl;
  Simplex_tree st5(std::move(st1));
  Simplex_tree st6(std::move(st2));

  // Cross check
  BOOST_CHECK(st5 == st6);
  BOOST_CHECK(st == st6);
  BOOST_CHECK(st5 == st);

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF MOVE ASSIGNMENT" << std::endl;

  // A swap is a copy ctor of a tmp value, then it uses move assignment
  std::swap(st3, st4);
  BOOST_CHECK(st3 == st4);
  BOOST_CHECK(st == st4);
  BOOST_CHECK(st3 == st);

}
