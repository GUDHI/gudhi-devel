/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/random_graph_generators.h>
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <array>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "random_graph_generators"
#include <boost/test/unit_test.hpp>

using Edge = std::array<int, 2>;
using Edge_range = std::vector<Edge>;

BOOST_AUTO_TEST_CASE( random_edges_limits ) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_edges_limits )" << std::endl;
  Edge_range edges = Gudhi::random_edges(0);
  std::cout << "Gudhi::random_edges(0).size() = " << edges.size() << std::endl;
  BOOST_CHECK(edges.size() == 0);

  edges = Gudhi::random_edges(1);
  std::cout << "Gudhi::random_edges(1).size() = " << edges.size() << std::endl;
  BOOST_CHECK(edges.size() == 0);

  edges = Gudhi::random_edges(2, 1.);
  std::cout << "Gudhi::random_edges(2).size() = " << edges.size() << std::endl;
  BOOST_CHECK(edges.size() == 1);

  edges = Gudhi::random_edges(15, 0.);
  std::cout << "Gudhi::random_edges(15, 0.).size() = " << edges.size() << std::endl;
  BOOST_CHECK(edges.size() == 0);

  edges = Gudhi::random_edges(15, 1.);
  std::cout << "Gudhi::random_edges(15, 1.).size() = " << edges.size() << std::endl;
  BOOST_CHECK(edges.size() == 105);
}

BOOST_AUTO_TEST_CASE( random_edges_mean ) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_edges_mean )" << std::endl;
  int nb_edges {0};
  for (int idx = 0; idx < 100; ++idx) {
    Edge_range edges = Gudhi::random_edges(15);
    nb_edges += edges.size();
#ifdef DEBUG_TRACES
    std::cout << edges.size() << ", ";
#endif  // DEBUG_TRACES
  }
#ifdef DEBUG_TRACES
  std::cout << "\n";
#endif  // DEBUG_TRACES
  std::cout << "Total number of edges for 100 x Gudhi::random_edges(15) [aka. 15% of 105] = " << nb_edges << std::endl;
  // 1575 +/- 10%
  BOOST_CHECK(nb_edges < 1733);
  BOOST_CHECK(nb_edges > 1417);
}

BOOST_AUTO_TEST_CASE( simplex_tree_random_flag_complex_test ) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( simplex_tree_random_flag_complex_test )" << std::endl;
  Gudhi::Simplex_tree<> stree;
  // Insert 100% of the possible edges, with 10 vertices
  simplex_tree_random_flag_complex(stree, 10, 1.);

  std::cout << "Random flag complex with " << stree.num_vertices() << " vertices and " << stree.num_simplices() <<
               " simplices" << std::endl;
  BOOST_CHECK(stree.num_vertices() == 10);
  BOOST_CHECK(stree.num_simplices() == 55);
}

BOOST_AUTO_TEST_CASE( simplex_tree_random_flag_complex_test_1000 ) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( simplex_tree_random_flag_complex_test_1000 )" << std::endl;
  Gudhi::Simplex_tree<> stree;
  // Insert 15% of the possible edges, with 1000 vertices
  simplex_tree_random_flag_complex(stree, 1000);

  std::cout << "Random flag complex with " << stree.num_vertices() << " vertices and " << stree.num_simplices() <<
               " simplices" << std::endl;

  // 15 % of 499500 = 74925
  // 74925 +/- 5% = [71178.75, 78671.25]
  BOOST_CHECK(stree.num_simplices() < 78671);
  BOOST_CHECK(stree.num_simplices() > 71178);
}
