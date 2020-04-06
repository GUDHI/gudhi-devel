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
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>  // std::pair, std::make_pair
#include <cmath>  // float comparison
#include <limits>
#include <functional>  // greater
#include <tuple>  // std::tie

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "collapse"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/FlagComplexSpMatrix.h"
#include "gudhi/Rips_edge_list.h"
#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/distance_functions.h"

//using namespace Gudhi;

// Types definition
//using Vector_of_points = std::vector<std::vector<double>>;

//using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
//using Filtration_value = double;
//using Rips_edge_list = Gudhi::rips_edge_list::Rips_edge_list<double>;
/*using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
*/
//using Distance_matrix = std::vector<std::vector<double>>;

using Filtration_value = double;
using Vertex_handle = size_t;
using Filtered_edge = std::tuple<Filtration_value, Vertex_handle, Vertex_handle>;
using Filtered_sorted_edge_list = std::vector<std::tuple<Filtration_value, Vertex_handle, Vertex_handle>>;

bool find_edge_in_list(const Filtered_edge& edge, const Filtered_sorted_edge_list& edge_list) {
  for (auto edge_from_list : edge_list) {
    if (edge_from_list == edge)
      return true;
  }
  return false;
}

void trace_and_check_collapse(const Filtered_sorted_edge_list& edges, const Filtered_sorted_edge_list& removed_edges) {
  std::cout << "BEFORE COLLAPSE - Total number of edges: " << edges.size() << std::endl;
  BOOST_CHECK(edges.size() > 0);
  for (auto edge : edges) {
    std::cout << "f[" << std::get<1>(edge) << ", " << std::get<2>(edge) << "] = " << std::get<0>(edge) << std::endl;
  }

  FlagComplexSpMatrix flag_complex_sparse_matrix(5, edges);
  auto collapse_edges = flag_complex_sparse_matrix.filtered_edge_collapse();
  std::cout << "AFTER COLLAPSE - Total number of edges: " << collapse_edges.size() << std::endl;
  BOOST_CHECK(collapse_edges.size() <= edges.size());
  for (auto edge_from_collapse : collapse_edges) {
    std::cout << "f[" << std::get<1>(edge_from_collapse) << ", " << std::get<2>(edge_from_collapse) << "] = "
              << std::get<0>(edge_from_collapse) << std::endl;
    // Check each edge from collapse is in the input
    BOOST_CHECK(find_edge_in_list(edge_from_collapse, edges));
  }

  std::cout << "CHECK COLLAPSE - Total number of removed edges: " << removed_edges.size() << std::endl;
  for (auto removed_edge : removed_edges) {
    std::cout << "f[" << std::get<1>(removed_edge) << ", " << std::get<2>(removed_edge) << "] = "
              << std::get<0>(removed_edge) << std::endl;
    // Check each removed edge from collapse is in the input
    BOOST_CHECK(!find_edge_in_list(removed_edge, collapse_edges));
  }

}

BOOST_AUTO_TEST_CASE(collapse) {
  /*
      1   2
      o---o
      |   |
      |   |
      |   |
      o---o
      0   3
  */
  Filtered_sorted_edge_list edges {{1., 0, 1}, {1., 1, 2}, {1., 2, 3}, {1., 3, 0}};
  trace_and_check_collapse(edges, {});
  
  /*
      1   2
      o---o
      |\ /|
      | x |
      |/ \|
      o---o
      0   3
  */
  edges.push_back({2., 0, 2});
  edges.push_back({2., 1, 3});
  trace_and_check_collapse(edges, {{2., 1, 3}});

  /*
      1   2   4
      o---o---o
      |\ /|   |
      | x |   |
      |/ \|   |
      o---o---o
      0   3   5
  */
  edges.push_back({3., 2, 4});
  edges.push_back({3., 4, 5});
  edges.push_back({3., 5, 3});
  trace_and_check_collapse(edges, {{2., 1, 3}});

  /*
      1   2   4
      o---o---o
      |\ /|\ /|
      | x | x |
      |/ \|/ \|
      o---o---o
      0   3   5
  */
  edges.push_back({4., 2, 5});
  edges.push_back({4., 4, 3});
  trace_and_check_collapse(edges, {{2., 1, 3}, {4., 2, 5}, {4., 4, 3}});

  /*
      1   2   4
      o---o---o
      |\ /|\ /|
      | x | x |  + [0,4] and [1,5]
      |/ \|/ \|
      o---o---o
      0   3   5
  */
  edges.push_back({5., 1, 5});
  edges.push_back({5., 0, 4});
  trace_and_check_collapse(edges, {{2., 1, 3}, {4., 2, 5}, {4., 4, 3}, {5., 0, 4}});
}


