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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_graph_expansion"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "gudhi/Simplex_tree.h"
#include <gudhi/Unitary_tests_utils.h>

using namespace Gudhi;

typedef boost::mpl::list<Simplex_tree<>, Simplex_tree<Simplex_tree_options_fast_persistence>> list_of_tested_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion_all_is_blocked, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_expansion_all_is_blocked\n";
  std::clog << "********************************************************************\n";
  using Simplex_handle = typename typeST::Simplex_handle;
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  typeST stree_copy = simplex_tree;

  simplex_tree.expansion_with_blockers(3, [&](Simplex_handle sh){ return true; });

  std::clog << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::clog << " - dimension " << simplex_tree.dimension() << "\n";
  std::clog << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::clog << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

  BOOST_CHECK(stree_copy == simplex_tree);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion_with_blockers_3, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_expansion_with_blockers_3\n";
  std::clog << "********************************************************************\n";
  using Simplex_handle = typename typeST::Simplex_handle;
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  simplex_tree.expansion_with_blockers(3, [&](Simplex_handle sh){
      bool result = false;
      std::clog << "Blocker on [";
      // User can loop on the vertices from the given simplex_handle i.e.
      for (auto vertex : simplex_tree.simplex_vertex_range(sh)) {
        // We block the expansion, if the vertex '6' is in the given list of vertices
        if (vertex == 6)
          result = true;
        std::clog << vertex << ", ";
      }
      std::clog << "] ( " << simplex_tree.filtration(sh);
      // User can re-assign a new filtration value directly in the blocker (default is the maximal value of boudaries)
      simplex_tree.assign_filtration(sh, simplex_tree.filtration(sh) + 1.);

      std::clog << " + 1. ) = " << result << std::endl;

      return result;
    });

  std::clog << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::clog << " - dimension " << simplex_tree.dimension() << "\n";
  std::clog << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::clog << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

  BOOST_CHECK(simplex_tree.num_simplices() == 23);
  BOOST_CHECK(simplex_tree.dimension() == 3);
  // {4, 5, 6} shall be blocked
  BOOST_CHECK(simplex_tree.find({4, 5, 6}) == simplex_tree.null_simplex());
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,2})),
                                                          static_cast<typename typeST::Filtration_value>(4.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,3})),
                                                          static_cast<typename typeST::Filtration_value>(5.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(6.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({1,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(6.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(7.));

}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion_with_blockers_2, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_expansion_with_blockers_2\n";
  std::clog << "********************************************************************\n";
  using Simplex_handle = typename typeST::Simplex_handle;
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  simplex_tree.expansion_with_blockers(2, [&](Simplex_handle sh){
      bool result = false;
      std::clog << "Blocker on [";
      // User can loop on the vertices from the given simplex_handle i.e.
      for (auto vertex : simplex_tree.simplex_vertex_range(sh)) {
        // We block the expansion, if the vertex '6' is in the given list of vertices
        if (vertex == 6)
          result = true;
        std::clog << vertex << ", ";
      }
      std::clog << "] ( " << simplex_tree.filtration(sh);
      // User can re-assign a new filtration value directly in the blocker (default is the maximal value of boudaries)
      simplex_tree.assign_filtration(sh, simplex_tree.filtration(sh) + 1.);

      std::clog << " + 1. ) = " << result << std::endl;

      return result;
    });

  std::clog << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::clog << " - dimension " << simplex_tree.dimension() << "\n";
  std::clog << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::clog << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

  BOOST_CHECK(simplex_tree.num_simplices() == 22);
  BOOST_CHECK(simplex_tree.dimension() == 2);
  // {4, 5, 6} shall be blocked
  BOOST_CHECK(simplex_tree.find({4, 5, 6}) == simplex_tree.null_simplex());
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,2})),
                                                          static_cast<typename typeST::Filtration_value>(4.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,3})),
                                                          static_cast<typename typeST::Filtration_value>(5.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(6.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({1,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(6.));
  BOOST_CHECK(simplex_tree.find({0,1,2,3}) == simplex_tree.null_simplex());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion_with_find_simplex_blockers, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_expansion_with_find_simplex_blockers\n";
  std::clog << "********************************************************************\n";
  using Simplex_handle = typename typeST::Simplex_handle;
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  simplex_tree.expansion_with_blockers(3, [&](Simplex_handle sh){
      bool result = false;
      std::clog << "Blocker on [";
      std::vector<typename typeST::Vertex_handle> simplex;
      // User can loop on the vertices from the given simplex_handle i.e.
      for (auto vertex : simplex_tree.simplex_vertex_range(sh)) {
        // We block the expansion, if the vertex '1' is in the given list of vertices
        if (vertex == 1)
          result = true;
        std::clog << vertex << ", ";
        simplex.push_back(vertex);
      }
      std::clog << "] => " << result << std::endl;
      // Not efficient but test it works - required by the python interface
      BOOST_CHECK(simplex_tree.find(simplex) == sh);
      return result;
    });

  std::clog << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::clog << " - dimension " << simplex_tree.dimension() << "\n";
  std::clog << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::clog << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

  BOOST_CHECK(simplex_tree.num_simplices() == 20);
  BOOST_CHECK(simplex_tree.dimension() == 2);

  // {1, 2, 3}, {0, 1, 2} and {0, 1, 3} shall be blocked as it contains vertex 1
  BOOST_CHECK(simplex_tree.find({4, 5, 6}) != simplex_tree.null_simplex());
  BOOST_CHECK(simplex_tree.find({1, 2, 3}) == simplex_tree.null_simplex());
  BOOST_CHECK(simplex_tree.find({0, 2, 3}) != simplex_tree.null_simplex());
  BOOST_CHECK(simplex_tree.find({0, 1, 2}) == simplex_tree.null_simplex());
  BOOST_CHECK(simplex_tree.find({0, 1, 3}) == simplex_tree.null_simplex());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion_3, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_expansion_3\n";
  std::clog << "********************************************************************\n";
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  simplex_tree.expansion(3);
  std::clog << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::clog << " - dimension " << simplex_tree.dimension() << "\n";
  std::clog << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::clog << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

  BOOST_CHECK(simplex_tree.num_simplices() == 24);
  BOOST_CHECK(simplex_tree.dimension() == 3);

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({4,5,6})),
                                                          static_cast<typename typeST::Filtration_value>(10.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,2})),
                                                          static_cast<typename typeST::Filtration_value>(3.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,3})),
                                                          static_cast<typename typeST::Filtration_value>(4.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(5.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({1,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(5.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(5.));

}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion_2, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_expansion_2\n";
  std::clog << "********************************************************************\n";
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  simplex_tree.expansion(2);

  std::clog << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::clog << " - dimension " << simplex_tree.dimension() << "\n";
  std::clog << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::clog << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

  BOOST_CHECK(simplex_tree.num_simplices() == 23);
  BOOST_CHECK(simplex_tree.dimension() == 2);

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({4,5,6})),
                                                          static_cast<typename typeST::Filtration_value>(10.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,2})),
                                                          static_cast<typename typeST::Filtration_value>(3.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,1,3})),
                                                          static_cast<typename typeST::Filtration_value>(4.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({0,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(5.));
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(simplex_tree.find({1,2,3})),
                                                          static_cast<typename typeST::Filtration_value>(5.));
  BOOST_CHECK(simplex_tree.find({0,1,2,3}) == simplex_tree.null_simplex());
}
