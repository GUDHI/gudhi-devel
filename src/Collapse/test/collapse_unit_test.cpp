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
#include <tuple>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "collapse"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "gudhi/Flag_complex_sparse_matrix.h"

using Filtration_value = float;
using Vertex_handle = short;
using Filtered_edge = std::tuple<Filtration_value, Vertex_handle, Vertex_handle>;
using Filtered_sorted_edge_list = std::vector<Filtered_edge>;
using Flag_complex_sparse_matrix = Gudhi::collapse::Flag_complex_sparse_matrix<Vertex_handle, Filtration_value>;

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

  std::cout << "COLLAPSE - keep edges: " << std::endl;
  Flag_complex_sparse_matrix flag_complex_sparse_matrix(edges);
  Filtered_sorted_edge_list collapse_edges;
  flag_complex_sparse_matrix.filtered_edge_collapse(
    [&collapse_edges](std::pair<Vertex_handle, Vertex_handle> edge, Filtration_value filtration) {
      std::cout << "f[" << std::get<0>(edge) << ", " << std::get<1>(edge) << "] = " << filtration << std::endl;
        collapse_edges.push_back({filtration, std::get<0>(edge), std::get<1>(edge)});
      });
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
  trace_and_check_collapse(edges, {{2., 1, 3}, {4., 4, 3}});

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
  trace_and_check_collapse(edges, {{2., 1, 3}, {4., 4, 3}, {5., 0, 4}});
}


