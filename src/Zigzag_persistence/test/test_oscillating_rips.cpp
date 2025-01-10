/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>
#include <algorithm>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "zigzag_persistence"
#include <boost/test/unit_test.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Zigzag_persistence/oscillating_rips_iterators.h>

struct Simplex_tree_options_oscillating_rips {
  typedef Gudhi::linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::int64_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = true;
  static const bool link_nodes_by_label = true;
  static const bool stable_simplex_handles = true;
};

using Gudhi::zigzag_persistence::Oscillating_rips_edge_order_policy;
using Gudhi::zigzag_persistence::Oscillating_rips_edge_range;
using Gudhi::zigzag_persistence::Oscillating_rips_simplex_range;
using Gudhi::zigzag_persistence::Zigzag_edge;

using Point = std::vector<double>;
using StableFilteredComplex = Gudhi::Simplex_tree<Simplex_tree_options_oscillating_rips>;
using Simplex_handle = typename StableFilteredComplex::Simplex_handle;
using Vertex_handle = typename StableFilteredComplex::Vertex_handle;
using Filtration_value = typename StableFilteredComplex::Filtration_value;
using Edge = Zigzag_edge<Filtration_value>;
using Face = std::vector<Vertex_handle>;
using Arrow = std::tuple<Simplex_handle, Filtration_value, bool>;
using VArrow = std::tuple<Face, Filtration_value, bool>;
using EdgeRange = Oscillating_rips_edge_range<Filtration_value>;
template <class EdgeRangeIterator>
using OscillatingRipsSimplexRange = Oscillating_rips_simplex_range<StableFilteredComplex, EdgeRangeIterator>;

std::vector<Point> get_point_cloud() { return {{0, 0}, {1.2, 1.2}, {0, 1.1}, {1, 0}, {0.51, 0.49}}; }

std::vector<double> get_epsilons(const std::vector<Point>& points)
{
  auto size = points.size();
  std::vector<Point> sortedPoints;
  sortedPoints.reserve(size);
  std::vector<Filtration_value> epsilonValues;
  epsilonValues.reserve(size);

  Gudhi::subsampling::choose_n_farthest_points(Gudhi::Euclidean_distance(),
                                               points,
                                               size,  // final size
                                               0,     // starting point
                                               std::back_inserter(sortedPoints),
                                               std::back_inserter(epsilonValues));
  // need to shift values output by subsampling:
  for (unsigned int i = 1; i < size; ++i) {
    epsilonValues[i - 1] = epsilonValues[i];
  }
  epsilonValues[size - 1] = 0;

  BOOST_CHECK(sortedPoints == points);

  return epsilonValues;
}

std::vector<VArrow> get_filtration(const std::vector<double>& eps)
{
  Filtration_value inf = std::numeric_limits<Filtration_value>::infinity();
  std::vector<VArrow> arrows;
  arrows.reserve(46);

  arrows.emplace_back(Face{0}, inf, true);

  arrows.emplace_back(Face{1}, eps[0], true);
  arrows.emplace_back(Face{0, 1}, eps[0], true);

  arrows.emplace_back(Face{2}, eps[1], true);
  arrows.emplace_back(Face{1, 2}, eps[1], true);
  arrows.emplace_back(Face{0, 2}, eps[1], true);
  arrows.emplace_back(Face{0, 1, 2}, eps[1], true);

  arrows.emplace_back(Face{3}, eps[2], true);
  arrows.emplace_back(Face{2, 3}, eps[2], true);
  arrows.emplace_back(Face{1, 3}, eps[2], true);
  arrows.emplace_back(Face{1, 2, 3}, eps[2], true);
  arrows.emplace_back(Face{0, 3}, eps[2], true);
  arrows.emplace_back(Face{0, 1, 3}, eps[2], true);
  arrows.emplace_back(Face{0, 2, 3}, eps[2], true);

  arrows.emplace_back(Face{0, 2, 3}, eps[2], false);
  arrows.emplace_back(Face{0, 1, 3}, eps[2], false);
  arrows.emplace_back(Face{1, 2, 3}, eps[2], false);
  arrows.emplace_back(Face{2, 3}, eps[2], false);
  arrows.emplace_back(Face{0, 1, 2}, eps[2], false);
  arrows.emplace_back(Face{0, 1}, eps[2], false);

  arrows.emplace_back(Face{4}, eps[3], true);
  arrows.emplace_back(Face{1, 4}, eps[3], true);
  arrows.emplace_back(Face{2, 4}, eps[3], true);
  arrows.emplace_back(Face{1, 2, 4}, eps[3], true);
  arrows.emplace_back(Face{0, 4}, eps[3], true);
  arrows.emplace_back(Face{0, 2, 4}, eps[3], true);
  arrows.emplace_back(Face{3, 4}, eps[3], true);
  arrows.emplace_back(Face{1, 3, 4}, eps[3], true);
  arrows.emplace_back(Face{0, 3, 4}, eps[3], true);

  arrows.emplace_back(Face{0, 3, 4}, eps[3], false);
  arrows.emplace_back(Face{1, 3, 4}, eps[3], false);
  arrows.emplace_back(Face{3, 4}, eps[3], false);
  arrows.emplace_back(Face{0, 2, 4}, eps[3], false);
  arrows.emplace_back(Face{0, 4}, eps[3], false);
  arrows.emplace_back(Face{1, 2, 4}, eps[3], false);
  arrows.emplace_back(Face{2, 4}, eps[3], false);
  arrows.emplace_back(Face{1, 4}, eps[3], false);
  arrows.emplace_back(Face{0, 3}, eps[3], false);
  arrows.emplace_back(Face{1, 3}, eps[3], false);
  arrows.emplace_back(Face{0, 2}, eps[3], false);
  arrows.emplace_back(Face{1, 2}, eps[3], false);

  arrows.emplace_back(Face{4}, -inf, false);
  arrows.emplace_back(Face{3}, -inf, false);
  arrows.emplace_back(Face{2}, -inf, false);
  arrows.emplace_back(Face{1}, -inf, false);
  arrows.emplace_back(Face{0}, -inf, false);

  return arrows;
}

template <class EdgeRangeIterator>
void test_edges(EdgeRangeIterator& start, const EdgeRangeIterator& end, const std::vector<double>& eps)
{
  auto comp = [](const Edge& e1, const Edge& e2) {
    if (e2.get_filtration_value() == e1.get_filtration_value()) {
      if (e1.get_direction() == e2.get_direction()) {
        if (e1.get_smallest_vertex() == e2.get_smallest_vertex()) {
          if (e1.get_biggest_vertex() == e2.get_biggest_vertex()) return false;
          return e1.get_biggest_vertex() < e2.get_biggest_vertex();
        } else {
          return e1.get_smallest_vertex() < e2.get_smallest_vertex();
        }
      } else {
        return e1.get_direction();
      }
    } else {
      return e1.get_filtration_value() > e2.get_filtration_value();
    }
  };

  std::vector<Edge> realEdges;
  realEdges.reserve(30);
  for (const VArrow& a : get_filtration(eps)) {
    const Face& f = std::get<0>(a);
    if (f.size() < 3) {
      if (f.size() == 1) {
        realEdges.emplace_back(f[0], f[0], std::get<1>(a), std::get<2>(a));
      } else if (f.size() == 2) {
        realEdges.emplace_back(f[0], f[1], std::get<1>(a), std::get<2>(a));
      }
    }
  }
  // for same filtration value, the order can change for the different iterator types.
  std::sort(realEdges.begin(), realEdges.end(), comp);

  std::vector<Edge> testEdges;
  testEdges.reserve(30);
  unsigned int i = 0;
  while (start != end && i < realEdges.size()) {  // to avoid infinite loop when something is wrong
    testEdges.push_back(*start);
    ++start;
    ++i;
  }
  BOOST_CHECK(start == end);
  // for same filtration value, the order can change for the different iterator types.
  std::sort(testEdges.begin(), testEdges.end(), comp);

  // instead of comparing directly testEdges == realEdges, because it makes debuging easier with more details
  i = 0;
  BOOST_CHECK(testEdges.size() == realEdges.size());
  for (const auto& e : testEdges) {
    // std::cout << "test: " << e << "\n";
    // std::cout << "real: " << realEdges[i] << "\n";
    BOOST_CHECK(e == realEdges[i]);
    ++i;
  }
}

template <class SimplexRangeIterator>
void test_filtration(SimplexRangeIterator& start,
                     const SimplexRangeIterator& end,
                     const StableFilteredComplex& st,
                     const std::vector<double>& eps)
{
  auto comp = [](const VArrow& e1, const VArrow& e2) {
    if (std::get<1>(e1) == std::get<1>(e2)) {
      if (std::get<2>(e1) == std::get<2>(e2)) {
        const auto& v1 = std::get<0>(e1);
        const auto& v2 = std::get<0>(e2);
        return std::lexicographical_compare(v1.begin(), v1.end(), v2.begin(), v2.end());
      } else {
        return std::get<2>(e1);
      }
    } else {
      return std::get<1>(e1) > std::get<1>(e2);
    }
  };

  std::vector<VArrow> realSimplices = get_filtration(eps);
  // for same filtration value, the order can change for the different iterator types.
  std::sort(realSimplices.begin(), realSimplices.end(), comp);

  std::vector<VArrow> testSimplices;
  testSimplices.reserve(46);
  unsigned int i = 0;
  while (start != end && i < realSimplices.size()) {  // to avoid infinite loop when something is wrong
    testSimplices.emplace_back();
    auto& vertices = std::get<0>(testSimplices.back());
    for (Vertex_handle v : st.simplex_vertex_range(std::get<0>(*start))) {
      vertices.push_back(v);
    }
    std::sort(vertices.begin(), vertices.end());
    std::get<1>(testSimplices.back()) = std::get<1>(*start);
    std::get<2>(testSimplices.back()) = std::get<2>(*start);
    ++start;
    ++i;
  }
  BOOST_CHECK(start == end);
  // for same filtration value, the order can change for the different iterator types.
  std::sort(testSimplices.begin(), testSimplices.end(), comp);

  i = 0;
  BOOST_CHECK(testSimplices.size() == realSimplices.size());
  for (const auto& e : testSimplices) {
    // std::cout << "test: ";
    // for (const auto& v : std::get<0>(e)) {
    //   std::cout << v << " ";
    // }
    // std::cout << "\n";
    // std::cout << "real: ";
    // for (const auto& v : std::get<0>(realSimplices[i])) {
    //   std::cout << v << " ";
    // }
    // std::cout << "\n";
    BOOST_CHECK(std::get<0>(e) == std::get<0>(realSimplices[i]));
    BOOST_CHECK(std::get<1>(e) == std::get<1>(realSimplices[i]));
    BOOST_CHECK(std::get<2>(e) == std::get<2>(realSimplices[i]));
    ++i;
  }
}

BOOST_AUTO_TEST_CASE(oscillating_rips_edge_range1)
{
  double nu = 1.73, mu = 2;
  std::vector<Point> points = get_point_cloud();
  std::vector<double> eps = {1.697, 1.1, 1, 0.707, 0};

  auto start1 = EdgeRange::begin(nu, mu, points, Gudhi::Euclidean_distance(), eps);
  test_edges(start1, EdgeRange::end(), eps);

  auto vec1 = EdgeRange::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), eps);
  auto start3 = vec1.begin();
  test_edges(start3, vec1.end(), eps);
}

BOOST_AUTO_TEST_CASE(oscillating_rips_edge_range2)
{
  double nu = 1.76, mu = 2;
  std::vector<Point> points = get_point_cloud();
  Oscillating_rips_edge_order_policy p = Oscillating_rips_edge_order_policy::FARTHEST_POINT_ORDERING;
  std::vector<double> eps = get_epsilons(points);

  auto start2 = EdgeRange::begin(nu, mu, points, Gudhi::Euclidean_distance(), p);
  test_edges(start2, EdgeRange::end(), eps);

  auto vec2 = EdgeRange::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);
  auto start4 = vec2.begin();
  test_edges(start4, vec2.end(), eps);
}

BOOST_AUTO_TEST_CASE(oscillating_rips_simplex_range1)
{
  double nu = 1.73, mu = 2;
  int maxDim = 2;
  std::vector<Point> points = get_point_cloud();
  std::vector<double> eps = {1.697, 1.1, 1, 0.707, 0};
  StableFilteredComplex st;

  auto edgeStart1 = EdgeRange::begin(nu, mu, points, Gudhi::Euclidean_distance(), eps);
  auto edgeEnd1 = EdgeRange::end();
  auto start1 = OscillatingRipsSimplexRange<typename EdgeRange::Oscillating_rips_edge_iterator>::begin(
      edgeStart1, edgeEnd1, st, maxDim);
  auto end1 = OscillatingRipsSimplexRange<typename EdgeRange::Oscillating_rips_edge_iterator>::end();
  test_filtration(start1, end1, st, eps);

  BOOST_CHECK(st.is_empty());

  auto vec1 = EdgeRange::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), eps);
  auto edgeStart3 = vec1.cbegin();
  auto edgeEnd3 = vec1.cend();
  auto start3 = OscillatingRipsSimplexRange<std::vector<Edge>::const_iterator>::begin(edgeStart3, edgeEnd3, st, maxDim);
  auto end3 = OscillatingRipsSimplexRange<std::vector<Edge>::const_iterator>::end();
  test_filtration(start3, end3, st, eps);
}

BOOST_AUTO_TEST_CASE(oscillating_rips_simplex_range2)
{
  double nu = 1.76, mu = 2;
  int maxDim = 2;
  std::vector<Point> points = get_point_cloud();
  Oscillating_rips_edge_order_policy p = Oscillating_rips_edge_order_policy::FARTHEST_POINT_ORDERING;
  std::vector<double> eps = get_epsilons(points);
  StableFilteredComplex st;

  auto edgeStart1 = EdgeRange::begin(nu, mu, points, Gudhi::Euclidean_distance(), p);
  auto edgeEnd1 = EdgeRange::end();
  auto start1 = OscillatingRipsSimplexRange<typename EdgeRange::Oscillating_rips_edge_iterator>::begin(
      edgeStart1, edgeEnd1, st, maxDim);
  auto end1 = OscillatingRipsSimplexRange<typename EdgeRange::Oscillating_rips_edge_iterator>::end();
  test_filtration(start1, end1, st, eps);

  BOOST_CHECK(st.is_empty());

  auto vec1 = EdgeRange::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);
  auto edgeStart3 = vec1.cbegin();
  auto edgeEnd3 = vec1.cend();
  auto start3 = OscillatingRipsSimplexRange<std::vector<Edge>::const_iterator>::begin(edgeStart3, edgeEnd3, st, maxDim);
  auto end3 = OscillatingRipsSimplexRange<std::vector<Edge>::const_iterator>::end();
  test_filtration(start3, end3, st, eps);
}
