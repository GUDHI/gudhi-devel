/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <limits>  // for NaN
#include <cmath>  // for isNaN

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_make_filtration_non_decreasing"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

typedef boost::mpl::list<Simplex_tree<>, Simplex_tree<Simplex_tree_options_fast_persistence>> list_of_tested_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(make_filtration_non_decreasing, typeST, list_of_tested_variants) {
  typeST st;

  st.insert_simplex_and_subfaces({2, 1, 0}, 2.0);
  st.insert_simplex_and_subfaces({3, 0}, 2.0);
  st.insert_simplex_and_subfaces({3, 4, 5}, 2.0);
  
  /* Inserted simplex:     */
  /*    1                  */
  /*    o                  */
  /*   /X\                 */
  /*  o---o---o---o        */
  /*  2   0   3\X/4        */
  /*            o          */
  /*            5          */

  std::clog << "Check default insertion ensures the filtration values are non decreasing" << std::endl;
  BOOST_CHECK(!st.make_filtration_non_decreasing());

  // Because of non decreasing property of simplex tree, { 0 } , { 1 } and { 0, 1 } are going to be set from value 2.0
  // to 1.0
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, 1.0);
  
  // Inserted simplex:
  //    1   6
  //    o---o
  //   /X\7/
  //  o---o---o---o
  //  2   0   3\X/4
  //            o
  //            5
  
  std::clog << "Check default second insertion ensures the filtration values are non decreasing" << std::endl;
  BOOST_CHECK(!st.make_filtration_non_decreasing());
  
  // Copy original simplex tree
  typeST st_copy = st;

  // Modify specific values for st to become like st_copy thanks to make_filtration_non_decreasing
  st.assign_filtration(st.find({0,1,6,7}), 0.8);
  st.assign_filtration(st.find({0,1,6}), 0.9);
  st.assign_filtration(st.find({0,6}), 0.6);
  st.assign_filtration(st.find({3,4,5}), 1.2);
  st.assign_filtration(st.find({3,4}), 1.1);
  st.assign_filtration(st.find({4,5}), 1.99);
  
  std::clog << "Check the simplex_tree is rolled back in case of decreasing filtration values" << std::endl;
  BOOST_CHECK(st.make_filtration_non_decreasing());
  BOOST_CHECK(st == st_copy);

  // Other simplex tree
  typeST st_other;
  st_other.insert_simplex_and_subfaces({2, 1, 0}, 3.0);  // This one is different from st
  st_other.insert_simplex_and_subfaces({3, 0}, 2.0);
  st_other.insert_simplex_and_subfaces({3, 4, 5}, 2.0);
  st_other.insert_simplex_and_subfaces({0, 1, 6, 7}, 1.0);

  // Modify specific values for st to become like st_other thanks to make_filtration_non_decreasing
  st.assign_filtration(st.find({2}), 3.0);
  // By modifying just the simplex {2}
  // {0,1,2}, {1,2} and {0,2} will be modified
  
  std::clog << "Check the simplex_tree is repaired in case of decreasing filtration values" << std::endl;
  BOOST_CHECK(st.make_filtration_non_decreasing());
  BOOST_CHECK(st == st_other);

  // Modify specific values for st still to be non-decreasing
  st.assign_filtration(st.find({0,1,2}), 10.0);
  st.assign_filtration(st.find({0,2}), 9.0);
  st.assign_filtration(st.find({0,1,6,7}), 50.0);
  st.assign_filtration(st.find({0,1,6}), 49.0);
  st.assign_filtration(st.find({0,1,7}), 48.0);
  // Other copy simplex tree
  typeST st_other_copy = st;
  
  std::clog << "Check the simplex_tree is not modified in case of non-decreasing filtration values" << std::endl;
  BOOST_CHECK(!st.make_filtration_non_decreasing());
  BOOST_CHECK(st == st_other_copy);
  
}

BOOST_AUTO_TEST_CASE_TEMPLATE(make_filtration_non_decreasing_on_nan_values, typeST, list_of_tested_variants) {
  typeST st;

  st.insert_simplex_and_subfaces({2, 1, 0}, std::numeric_limits<double>::quiet_NaN());
  st.insert_simplex_and_subfaces({3, 0},    std::numeric_limits<double>::quiet_NaN());
  st.insert_simplex_and_subfaces({3, 4, 5}, std::numeric_limits<double>::quiet_NaN());
  
  /* Inserted simplex:     */
  /*    1                  */
  /*    o                  */
  /*   /X\                 */
  /*  o---o---o---o        */
  /*  2   0   3\X/4        */
  /*            o          */
  /*            5          */

  std::clog << "SPECIFIC CASE:" << std::endl;
  std::clog << "Insertion with NaN values does not ensure the filtration values are non decreasing" << std::endl;
  st.make_filtration_non_decreasing();

  std::clog << "Check all filtration values are NaN" << std::endl;
  for (auto f_simplex : st.complex_simplex_range()) {
    BOOST_CHECK(std::isnan(st.filtration(f_simplex)));
  }

  st.assign_filtration(st.find({0}), 0.);
  st.assign_filtration(st.find({1}), 0.);
  st.assign_filtration(st.find({2}), 0.);
  st.assign_filtration(st.find({3}), 0.);
  st.assign_filtration(st.find({4}), 0.);
  st.assign_filtration(st.find({5}), 0.);

  std::clog << "Check make_filtration_non_decreasing is modifying the simplicial complex" << std::endl;
  BOOST_CHECK(st.make_filtration_non_decreasing());
  
  std::clog << "Check all filtration values are now defined" << std::endl;
  for (auto f_simplex : st.complex_simplex_range()) {
    BOOST_CHECK(!std::isnan(st.filtration(f_simplex)));
  }
}
