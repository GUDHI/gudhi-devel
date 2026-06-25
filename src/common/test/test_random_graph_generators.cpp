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
#include <cmath>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "random_graph_generators"
#include <boost/test/unit_test.hpp>

using Edge = std::array<int, 2>;
using Edge_range = std::vector<Edge>;

BOOST_AUTO_TEST_CASE( random_edges_limits ) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_edges_limits )\n";
  Edge_range edges = Gudhi::random::random_edges(0);
  std::cout << "Gudhi::random_edges(0).size() = " << edges.size() << "\n";
  BOOST_CHECK(edges.size() == 0);

  edges = Gudhi::random::random_edges(1);
  std::cout << "Gudhi::random_edges(1).size() = " << edges.size() << "\n";
  BOOST_CHECK(edges.size() == 0);

  edges = Gudhi::random::random_edges(2, 1.);
  std::cout << "Gudhi::random_edges(2).size() = " << edges.size() << "\n";
  BOOST_CHECK(edges.size() == 1);

  edges = Gudhi::random::random_edges(15, 0.);
  std::cout << "Gudhi::random_edges(15, 0.).size() = " << edges.size() << "\n";
  BOOST_CHECK(edges.size() == 0);

  edges = Gudhi::random::random_edges(15, 1.);
  std::cout << "Gudhi::random_edges(15, 1.).size() = " << edges.size() << "\n";
  // when density is 1., it returns all possible edges. nb_edges = (nb_vertices * (nb_vertices - 1)) / 2
  BOOST_CHECK(edges.size() == (15 * 14) / 2);
  
  // Check edges vertex_handle
  for (auto edge: edges) {
    BOOST_CHECK(edge[0] < 15);
    BOOST_CHECK(edge[1] < 15);
  }
}

BOOST_AUTO_TEST_CASE( random_edges_density ) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( random_edges_density )\n";
  // default density is 0.15, it returns 15% of nb_edges_max = (nb_vertices * (nb_vertices - 1)) / 2
  Edge_range edges = Gudhi::random::random_edges(15);
  std::size_t nb_edges = std::round(0.15 * ((15 * 14) / 2));
  std::cout << "Total number of edges for Gudhi::random_edges(15) [aka. 15% of 105] = " << nb_edges << "\n";
  std::cout << "Gudhi::random_edges(15) returns " << edges.size() << " edges.\n";
  BOOST_CHECK(nb_edges == edges.size());

  // with density = 0.2, it returns 20% of nb_edges_max = (nb_vertices * (nb_vertices - 1)) / 2
  edges = Gudhi::random::random_edges(15, 0.2);
  nb_edges = std::round(0.2 * ((15 * 14) / 2));
  std::cout << "Total number of edges for Gudhi::random_edges(15, 0.2) [aka. 20% of 105] = " << nb_edges << "\n";
  std::cout << "Gudhi::random_edges(15, 0.2) returns " << edges.size() << " edges.\n";
}

BOOST_AUTO_TEST_CASE( simplex_tree_random_graph_test ) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( simplex_tree_random_graph_test )\n";
  Gudhi::Simplex_tree<> stree;
  // Insert 100% of the possible edges, with 10 vertices
  Gudhi::random::simplex_tree_random_graph(stree, 10, 1.);

  std::cout << "Random graph with " << stree.num_vertices() << " vertices and " << stree.num_simplices() <<
               " simplices\n";
  BOOST_CHECK(stree.num_vertices() == 10);
  BOOST_CHECK(stree.num_simplices() == 55);
}

BOOST_AUTO_TEST_CASE( simplex_tree_random_graph_test_1000 ) {
  std::cout << "  ## BOOST_AUTO_TEST_CASE( simplex_tree_random_graph_test_1000 )\n";
  Gudhi::Simplex_tree<> stree;
  // Insert 15% of the possible edges, with 1000 vertices
  Gudhi::random::simplex_tree_random_graph(stree, 1000);

  std::cout << "Random graph with " << stree.num_vertices() << " vertices and " << stree.num_simplices() <<
               " simplices\n";

  std::size_t expected_nb_edges = std::round(0.15 * ((1000*999)/2));
  std::cout << "(number of simplices - number of vertices) should be (0.15 * ((1000*999)/2) = "
            << expected_nb_edges << "\n";
  BOOST_CHECK((stree.num_simplices() - stree.num_vertices()) == expected_nb_edges);
  
  for (auto spx : stree.filtration_simplex_range()) {
    BOOST_CHECK(stree.filtration(spx) <= 1.);
    BOOST_CHECK(stree.filtration(spx) >= 0.);
  }

  // Empty the Simplex_tree
  stree.prune_above_dimension(-1);
  // Insert 20% of the possible edges, with 1000 vertices, and filtration values all equals to 0.
  Gudhi::random::simplex_tree_random_graph(stree, 1000, 0.2, 0., 0.);
  expected_nb_edges = std::round(0.2 * ((1000*999)/2));
  std::cout << "(number of simplices - number of vertices) should be (0.2 * ((1000*999)/2) = "
            << expected_nb_edges << "\n";
  BOOST_CHECK((stree.num_simplices() - stree.num_vertices()) == expected_nb_edges);
  for (auto spx : stree.filtration_simplex_range()) {
    BOOST_CHECK(stree.filtration(spx) == 0.);
  }

}
