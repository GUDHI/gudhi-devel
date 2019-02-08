/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2019 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "flag_complex_sparse_matrix"
#include <boost/test/unit_test.hpp>

#include <vector>
#include <tuple>

#include <gudhi/Flag_complex_sparse_matrix.h>
#include <gudhi/Unitary_tests_utils.h>

// Type definitions
using Filtered_edge = std::tuple<double, int, int>;
using Filtered_sorted_edge_list = std::vector<Filtered_edge>;
using Edge_list = std::vector<std::pair<int, int>>;
using Reduction_map = std::unordered_map<int, int>;

// Can be used with Edge_list or Reduction_map
template<typename Structure >
bool find(const Structure& structure, int edge_source, int edge_target) {
  bool found = false;
  for (auto structure_it:structure) {
    if ((structure_it.first == edge_source && structure_it.second == edge_target) ||
        (structure_it.first == edge_target && structure_it.second == edge_source))
      return true;
  }
  return found;
}

BOOST_AUTO_TEST_CASE(flag_complex_sparse_matrix_strong_collapse) {
  Filtered_sorted_edge_list input_edges;
  input_edges.push_back({1., 1, 2});
  input_edges.push_back({1., 2, 3});
  input_edges.push_back({1., 3, 4});
  input_edges.push_back({1., 4, 1});

  {
    // Edges:
    //  1 o---o 2
    //    |   |
    //    |   |
    //    |   |
    //  4 o---o 3
    // Check no collapse is performed
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(4, input_edges);

    mat_coll.strong_collapse();
    Edge_list collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
    BOOST_CHECK(collapsed_edges.size() == 8);
    BOOST_CHECK(find(collapsed_edges, 1, 2));
    BOOST_CHECK(find(collapsed_edges, 2, 3));
    BOOST_CHECK(find(collapsed_edges, 3, 4));
    BOOST_CHECK(find(collapsed_edges, 4, 1));

    // Added by strong_collapse method
    BOOST_CHECK(find(collapsed_edges, 1, 1));
    BOOST_CHECK(find(collapsed_edges, 2, 2));
    BOOST_CHECK(find(collapsed_edges, 3, 3));
    BOOST_CHECK(find(collapsed_edges, 4, 4));

    Reduction_map reduction_map = mat_coll.reduction_map();
    for (auto edge : reduction_map)
      std::cout << edge.first << " -> " << edge.second << std::endl;

    // No contraction is possible
    BOOST_CHECK(reduction_map.empty());

  }
  input_edges.push_back({1.5, 1, 3});
  input_edges.push_back({1.5, 2, 4});
  {
    // Edges:
    //  1 o---o 2
    //    |\ /|
    //    | X |
    //    |/ \|
    //  4 o---o 3
    // Check all is collapsed on 1
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(4, input_edges);

    mat_coll.strong_collapse();
    Edge_list collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;

    BOOST_CHECK(collapsed_edges.size() == 1);

    // Added by strong_collapse method
    BOOST_CHECK(find(collapsed_edges, 1, 1));

    Reduction_map reduction_map = mat_coll.reduction_map();
    for (auto edge : reduction_map)
      std::cout << edge.first << " -> " << edge.second << std::endl;

    BOOST_CHECK(find(reduction_map, 2, 1));
    BOOST_CHECK(find(reduction_map, 3, 1));
    BOOST_CHECK(find(reduction_map, 4, 1));
  }
  input_edges.push_back({2., 2, 5});
  input_edges.push_back({2., 5, 6});
  input_edges.push_back({2., 3, 6});
  {
    // Edges:
    //        2
    //  1 o---o---o 5
    //    |\ /|   |
    //    | X |   |
    //    |/ \|   |
    //  4 o---o---o 6
    //        3
    // Check 4 and 1 are collapsed on 2
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(6, input_edges);

    mat_coll.strong_collapse();
    Edge_list collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
    BOOST_CHECK(find(collapsed_edges, 2, 3));
    BOOST_CHECK(find(collapsed_edges, 2, 5));
    BOOST_CHECK(find(collapsed_edges, 3, 6));
    BOOST_CHECK(find(collapsed_edges, 5, 6));

    // Added by strong_collapse method
    BOOST_CHECK(find(collapsed_edges, 2, 2));
    BOOST_CHECK(find(collapsed_edges, 3, 3));
    BOOST_CHECK(find(collapsed_edges, 5, 5));
    BOOST_CHECK(find(collapsed_edges, 6, 6));

    Reduction_map reduction_map = mat_coll.reduction_map();
    for (auto edge : reduction_map)
      std::cout << edge.first << " -> " << edge.second << std::endl;

    BOOST_CHECK(find(reduction_map, 1, 2));
    BOOST_CHECK(find(reduction_map, 4, 2));
  }
  input_edges.push_back({2.5, 2, 6});
  input_edges.push_back({2.5, 5, 3});
  {
    // Edges:
    //        2
    //  1 o---o---o 5
    //    |\ /|\ /|
    //    | X | X |
    //    |/ \|/ \|
    //  4 o---o---o 6
    //        3
    // Check all is collapsed on 2
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(6, input_edges);

    mat_coll.strong_collapse();
    Edge_list collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
    // Added by strong_collapse method
    BOOST_CHECK(find(collapsed_edges, 2, 2));

    Reduction_map reduction_map = mat_coll.reduction_map();
    for (auto edge : reduction_map)
      std::cout << edge.first << " -> " << edge.second << std::endl;
    BOOST_CHECK(find(reduction_map, 1, 2));
    BOOST_CHECK(find(reduction_map, 3, 2));
    BOOST_CHECK(find(reduction_map, 4, 2));
    BOOST_CHECK(find(reduction_map, 5, 2));
    BOOST_CHECK(find(reduction_map, 6, 2));
  }
  input_edges.push_back({3., 3, 7});
  input_edges.push_back({3., 7, 8});
  input_edges.push_back({3., 4, 8});
  {
    // Edges:
    //        2
    //  1 o---o---o 5
    //    |\ /|\ /|
    //    | X | X |
    //    |/ \|/ \|
    //  4 o---o---o 6
    //    |   | 3
    //    |   |
    //    |   |
    //  8 o---o 7
    // Check 1, 2, 5 and 6 are collapsed on 3
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(8, input_edges);

    mat_coll.strong_collapse();
    Edge_list collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;

    BOOST_CHECK(find(collapsed_edges, 3, 4));
    BOOST_CHECK(find(collapsed_edges, 3, 7));
    BOOST_CHECK(find(collapsed_edges, 4, 8));
    BOOST_CHECK(find(collapsed_edges, 7, 8));

    // Added by strong_collapse method
    BOOST_CHECK(find(collapsed_edges, 3, 3));
    BOOST_CHECK(find(collapsed_edges, 4, 4));
    BOOST_CHECK(find(collapsed_edges, 7, 7));
    BOOST_CHECK(find(collapsed_edges, 8, 8));

    Reduction_map reduction_map = mat_coll.reduction_map();
    for (auto edge : reduction_map)
      std::cout << edge.first << " -> " << edge.second << std::endl;
    BOOST_CHECK(find(reduction_map, 1, 3));
    BOOST_CHECK(find(reduction_map, 2, 3));
    BOOST_CHECK(find(reduction_map, 5, 3));
    BOOST_CHECK(find(reduction_map, 6, 3));
  }
  input_edges.push_back({3., 3, 8});
  input_edges.push_back({3., 4, 7});
  {
    // Edges:
    //        2
    //  1 o---o---o 5
    //    |\ /|\ /|
    //    | X | X |
    //    |/ \|/ \|
    //  4 o---o---o 6
    //    |\ /| 3
    //    | X |
    //    |/ \|
    //  8 o---o 7
    // Check all is collapsed on 3
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(8, input_edges);

    mat_coll.strong_collapse();
    Edge_list collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;

    // Added by strong_collapse method
    BOOST_CHECK(find(collapsed_edges, 3, 3));

    Reduction_map reduction_map = mat_coll.reduction_map();
    for (auto edge : reduction_map)
      std::cout << edge.first << " -> " << edge.second << std::endl;
    BOOST_CHECK(find(reduction_map, 1, 3));
    BOOST_CHECK(find(reduction_map, 2, 3));
    BOOST_CHECK(find(reduction_map, 4, 3));
    BOOST_CHECK(find(reduction_map, 5, 3));
    BOOST_CHECK(find(reduction_map, 6, 3));
    BOOST_CHECK(find(reduction_map, 7, 3));
    BOOST_CHECK(find(reduction_map, 8, 3));
  }
}

BOOST_AUTO_TEST_CASE(flag_complex_sparse_matrix_dummy_test) {
  Filtered_sorted_edge_list input_edges;
  input_edges.push_back({1., 0, 1});
  input_edges.push_back({1., 1, 0});

  {
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(2, input_edges);
    mat_coll.print_sparse_skeleton();
    Edge_list collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
    Reduction_map reduction_map = mat_coll.reduction_map();
    for (auto edge : reduction_map)
      std::cout << edge.first << " -> " << edge.second << std::endl;

  }
}