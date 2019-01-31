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

BOOST_AUTO_TEST_CASE(flag_complex_sparse_matrix_strong_collapse) {
  Filtered_sorted_edge_list input_edges;
  input_edges.push_back({1., 1, 2});
  input_edges.push_back({1., 2, 3});
  input_edges.push_back({1., 3, 4});
  input_edges.push_back({1., 4, 1});

  {
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(4, input_edges);

    mat_coll.strong_collapse();
    auto collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
  }
  input_edges.push_back({1.5, 1, 3});
  input_edges.push_back({1.5, 2, 4});
  {
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(4, input_edges);

    mat_coll.strong_collapse();
    auto collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
  }
  input_edges.push_back({2., 2, 5});
  input_edges.push_back({2., 5, 6});
  input_edges.push_back({2., 3, 6});
  {
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(6, input_edges);

    mat_coll.strong_collapse();
    auto collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
  }
  input_edges.push_back({2.5, 2, 6});
  input_edges.push_back({2.5, 5, 3});
  {
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(6, input_edges);

    mat_coll.strong_collapse();
    auto collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
  }
  input_edges.push_back({3., 3, 7});
  input_edges.push_back({3., 7, 8});
  input_edges.push_back({3., 4, 8});
  {
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(8, input_edges);

    mat_coll.strong_collapse();
    auto collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
  }
  input_edges.push_back({3.5, 4, 7});
  input_edges.push_back({3.5, 3, 8});
  {
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(8, input_edges);

    mat_coll.strong_collapse();
    auto collapsed_edges = mat_coll.all_edges();
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    for (auto edge : collapsed_edges)
      std::cout << edge.first << " - " << edge.second << std::endl;
  }

}

