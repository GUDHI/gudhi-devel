/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "collapse"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <gudhi/Flag_complex_edge_collapser.h>
#include <gudhi/distance_functions.h>
#include <gudhi/graph_simplicial_complex.h>

#include <iostream>
#include <tuple>
#include <vector>
#include <array>
#include <cmath>

struct Simplicial_complex {
  using Vertex_handle = short;
  using Filtration_value = float;
};

using Vertex_handle = Simplicial_complex::Vertex_handle;
using Filtration_value = Simplicial_complex::Filtration_value;
using Filtered_edge = std::tuple<Vertex_handle, Vertex_handle, Filtration_value>;
using Filtered_edge_list = std::vector<Filtered_edge>;

template<typename Filtered_edge_range>
bool find_edge_in_list(const Filtered_edge& edge, const Filtered_edge_range& edge_list) {
  for (auto edge_from_list : edge_list) {
    if (edge_from_list == edge)
      return true;
  }
  return false;
}

template<typename Filtered_edge_range>
void trace_and_check_collapse(const Filtered_edge_range& filtered_edges, const Filtered_edge_list& removed_edges) {
  std::cout << "BEFORE COLLAPSE - Total number of edges: " << filtered_edges.size() << std::endl;
  BOOST_CHECK(filtered_edges.size() > 0);
  for (auto filtered_edge : filtered_edges) {
    std::cout << "f[" << std::get<0>(filtered_edge) << ", " << std::get<1>(filtered_edge) << "] = "
              << std::get<2>(filtered_edge) << std::endl;
  }

  std::cout << "COLLAPSE - keep edges: " << std::endl;
  auto remaining_edges = Gudhi::collapse::flag_complex_collapse_edges(filtered_edges);

  std::cout << "AFTER COLLAPSE - Total number of edges: " << remaining_edges.size() << std::endl;
  BOOST_CHECK(remaining_edges.size() <= filtered_edges.size());
  for (auto filtered_edge_from_collapse : remaining_edges) {
    std::cout << "f[" << std::get<0>(filtered_edge_from_collapse) << ", " << std::get<1>(filtered_edge_from_collapse)
              << "] = " << std::get<2>(filtered_edge_from_collapse) << std::endl;
    // Check each edge from collapse is in the input
    BOOST_CHECK(find_edge_in_list(filtered_edge_from_collapse, filtered_edges));
  }

  std::cout << "CHECK COLLAPSE - Total number of removed edges: " << removed_edges.size() << std::endl;
  for (auto removed_filtered_edge : removed_edges) {
    std::cout << "f[" << std::get<0>(removed_filtered_edge) << ", " << std::get<1>(removed_filtered_edge) << "] = "
              << std::get<2>(removed_filtered_edge) << std::endl;
    // Check each removed edge from collapse is in the input
    BOOST_CHECK(!find_edge_in_list(removed_filtered_edge, remaining_edges));
  }

}

BOOST_AUTO_TEST_CASE(collapse) {
  std::cout << "***** COLLAPSE *****" << std::endl;
  //  1   2
  //  o---o
  //  |   |
  //  |   |
  //  |   |
  //  o---o
  //  0   3
  Filtered_edge_list edges {{0, 1, 1.},
                            {1, 2, 1.},
                            {2, 3, 1.},
                            {3, 0, 1.}};
  trace_and_check_collapse(edges, {});
  
  //  1   2
  //  o---o
  //  |\ /|
  //  | x |
  //  |/ \|
  //  o---o
  //  0   3
  edges.emplace_back(0, 2, 2.);
  edges.emplace_back(1, 3, 2.);
  trace_and_check_collapse(edges, {{1, 3, 2.}});

  //  1   2   4
  //  o---o---o
  //  |\ /|   |
  //  | x |   |
  //  |/ \|   |
  //  o---o---o
  //  0   3   5
  edges.emplace_back(2, 4, 3.);
  edges.emplace_back(4, 5, 3.);
  edges.emplace_back(5, 3, 3.);
  trace_and_check_collapse(edges, {{1, 3, 2.}});

  //  1   2   4
  //  o---o---o
  //  |\ /|\ /|
  //  | x | x |
  //  |/ \|/ \|
  //  o---o---o
  //  0   3   5
  edges.emplace_back(2, 5, 4.);
  edges.emplace_back(4, 3, 4.);
  trace_and_check_collapse(edges, {{1, 3, 2.}, {4, 3, 4.}});

  //  1   2   4
  //  o---o---o
  //  |\ /|\ /|
  //  | x | x |  + [0,4] and [1,5]
  //  |/ \|/ \|
  //  o---o---o
  //  0   3   5
  edges.emplace_back(1, 5, 5.);
  edges.emplace_back(0, 4, 5.);
  trace_and_check_collapse(edges, {{1, 3, 2.}, {4, 3, 4.}, {0, 4, 5.}});
}

BOOST_AUTO_TEST_CASE(collapse_from_array) {
  std::cout << "***** COLLAPSE FROM ARRAY *****" << std::endl;
  //  1   2
  //  o---o
  //  |\ /|
  //  | x |
  //  |/ \|
  //  o---o
  //  0   3
  std::array<Filtered_edge, 6> f_edge_array = {{{0, 1, 1.},
                                                {1, 2, 1.},
                                                {2, 3, 1.},
                                                {3, 0, 1.},
                                                {0, 2, 2.},
                                                {1, 3, 2.}}};
  trace_and_check_collapse(f_edge_array, {{1, 3, 2.}});
}

BOOST_AUTO_TEST_CASE(collapse_from_proximity_graph) {
  std::cout << "***** COLLAPSE FROM PROXIMITY GRAPH *****" << std::endl;
  //  1   2
  //  o---o
  //  |\ /|
  //  | x |
  //  |/ \|
  //  o---o
  //  0   3
  std::vector<std::vector<Filtration_value>> point_cloud = {{0., 0.},
                                                            {0., 1.},
                                                            {1., 0.},
                                                            {1., 1.} };

  Filtration_value threshold = std::numeric_limits<Filtration_value>::infinity();
  using Proximity_graph = Gudhi::Proximity_graph<Simplicial_complex>;
  Proximity_graph proximity_graph = Gudhi::compute_proximity_graph<Simplicial_complex>(point_cloud,
                                                                                       threshold,
                                                                                       Gudhi::Euclidean_distance());

  auto remaining_edges = Gudhi::collapse::flag_complex_collapse_edges(
      boost::adaptors::transform(edges(proximity_graph), [&](auto&&edge){
        return std::make_tuple(static_cast<Vertex_handle>(source(edge, proximity_graph)),
                               static_cast<Vertex_handle>(target(edge, proximity_graph)),
                               get(Gudhi::edge_filtration_t(), proximity_graph, edge));
      })
    );

  BOOST_CHECK(remaining_edges.size() == 5);

  std::size_t filtration_is_edge_length_nb = 0;
  std::size_t filtration_is_diagonal_length_nb = 0;
  float epsilon = std::numeric_limits<Filtration_value>::epsilon();
  for (auto filtered_edge : remaining_edges) {
    if (std::get<2>(filtered_edge) == 1.)
      filtration_is_edge_length_nb++;
    if (std::fabs(std::get<2>(filtered_edge) - std::sqrt(2.)) <= epsilon)
      filtration_is_diagonal_length_nb++;
  }
  BOOST_CHECK(filtration_is_edge_length_nb == 4);
  BOOST_CHECK(filtration_is_diagonal_length_nb == 1);
}
