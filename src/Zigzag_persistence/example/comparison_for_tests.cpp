/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>  // for pair
#include <vector>
#include <unordered_set>

#include <boost/range/adaptors.hpp>

#include <gudhi/Zigzag_persistence/oscillating_rips_iterators.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Clock.h>

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_oscillating_rips>;
using Filtration_value = ST::Filtration_value;
using Square = Gudhi::zigzag_persistence::Square_root_edge_modifier<Filtration_value>;
using ZE = Gudhi::zigzag_persistence::Zigzag_edge<Filtration_value>;
using ORE = Gudhi::zigzag_persistence::Oscillating_rips_edge_range<Filtration_value>;
using ORi = Gudhi::zigzag_persistence::Oscillating_rips_simplex_range<ST,ORE::Oscillating_rips_edge_iterator>;
using ORv = Gudhi::zigzag_persistence::Oscillating_rips_simplex_range<ST,std::vector<ZE>::iterator>;
using Point = std::vector<double>;

void print_points(const std::vector<Point>& points) {
  std::cout << "Number of points: " << points.size() << "\n";
  for (const Point& p : points) {
    std::cout << "(" << p[0] << ", " << p[1] << ")\n";
  }
  std::cout << "\n";
}

std::vector<Point> build_point_cloud(unsigned int numberOfPoints, int seed) {
  std::vector<Point> finalPoints;
  std::set<Point> points;
  std::random_device dev;
  std::mt19937 rng(dev());
  if (seed > -1) rng.seed(seed);
  std::uniform_real_distribution<double> dist(0, 10);

  for (unsigned int i = 0; i < numberOfPoints; ++i) {
    auto res = points.insert({dist(rng), dist(rng)});
    while (!res.second) {
      res = points.insert({dist(rng), dist(rng)});
    }
    finalPoints.push_back(*res.first);
  }

  // print_points(finalPoints);

  return finalPoints;
}

void canonical_sort_edge(std::vector<ZE>& edges) {
  // canonical sort of the edges: as much as possible, edges should be removed in
  // the reverse order of their insertion. We decide to insert shorted edges first,
  // with increasing lexicographical order, and remove larger edges first, with
  // decreasing lexicographic order.

  // filtration then dimension, then lex order for insertion
  auto edge_cmp = [](const ZE& e1, const ZE& e2) {
    if (e1.get_filtration_value() != e2.get_filtration_value()) {
      return e1.get_filtration_value() < e2.get_filtration_value();
    }  // lower fil first

    if (e1.get_smallest_vertex() == e1.get_biggest_vertex()) {  // e1 is a vertex, -> put vertices first
      if (e2.get_smallest_vertex() == e2.get_biggest_vertex()) {
        return e1.get_smallest_vertex() < e2.get_smallest_vertex();
      }  //-> vertex of lower label
      else {
        return true;
      }       //-> always vertices before edges
    } else {  // e1 is an edge
      if (e2.get_smallest_vertex() == e2.get_biggest_vertex()) {
        return false;
      }       // e2 vertex, -> put it first
      else {  // both are edges, lexigraphic compare
        if (e1.get_smallest_vertex() != e2.get_smallest_vertex()) {
          return e1.get_smallest_vertex() < e2.get_smallest_vertex();
        }  // lex order
        if (e1.get_biggest_vertex() != e2.get_biggest_vertex()) {
          return e1.get_biggest_vertex() < e2.get_biggest_vertex();
        }
        return false;  // equality
      }
    }
  };
  // the inverse ordering for deletions
  auto inv_edge_cmp = [&](const ZE& e1, const ZE& e2) {
    if (e1.get_smallest_vertex() == e2.get_smallest_vertex() && e1.get_biggest_vertex() == e2.get_biggest_vertex()) {
      return false;
    }                            //== => false
    return !(edge_cmp(e1, e2));  // reverse order
  };
  // sort sequences of inclusions of same filtration with edge_cmp
  // sort sequences of removals of same filtration with inv_edge_cmp
  auto beg = edges.begin();
  auto end = edges.begin();
  auto curr_fil = beg->get_filtration_value();
  auto curr_type = beg->get_direction();
  while (beg != edges.end()) {
    while (end != edges.end() && end->get_filtration_value() == curr_fil && end->get_direction() == curr_type) {
      ++end;
    }
    if (curr_type) {
      sort(beg, end, edge_cmp);
    }  // sequence of insertions
    else {
      sort(beg, end, inv_edge_cmp);
    }  // sequence of removals
    beg = end;
    curr_fil = beg->get_filtration_value();
    curr_type = beg->get_direction();
  }
}

void test_edges_comp(const std::vector<Point>& points, double nu, double mu, ORE::Order_policy p) 
{
  std::vector<ZE> edges_v1 = ORE::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);

  unsigned int i = 0;
  for (const auto& e : ORE::get_iterator_range(nu, mu, points, Gudhi::Euclidean_distance(), p)) {
    if (i < edges_v1.size()) {
      if (!(edges_v1[i] == e)) {
        std::cout << "[" << i << "] different:\n";
        std::cout << edges_v1[i].get_smallest_vertex() << ", " << edges_v1[i].get_biggest_vertex() << ", "
                  << edges_v1[i].get_filtration_value() << ", " << edges_v1[i].get_direction() << "\n";
        std::cout << e.get_smallest_vertex() << ", " << e.get_biggest_vertex() << ", " << e.get_filtration_value()
                  << ", " << e.get_direction() << "\n";
      }/*  else {
        std::cout << "[" << i << "] same:\n";
        std::cout << edges_v1[i].get_smallest_vertex() << ", " << edges_v1[i].get_biggest_vertex() << ", "
                  << edges_v1[i].get_filtration_value() << ", " << edges_v1[i].get_direction() << "\n";
        std::cout << e.get_smallest_vertex() << ", " << e.get_biggest_vertex() << ", " << e.get_filtration_value()
                  << ", " << e.get_direction() << "\n";
      } */
    } else {
      std::cout << "[" << i << "] too long:\n";
      std::cout << e.get_smallest_vertex() << ", " << e.get_biggest_vertex() << ", " << e.get_filtration_value() << ", "
                << e.get_direction() << "\n";
    }
    ++i;
  }
}

void test_edges_canonical_sort(const std::vector<Point>& points, double nu, double mu, ORE::Order_policy p) {
  std::vector<ZE> edges_v1 = ORE::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);
  std::vector<ZE> ordered_edges_v1(edges_v1);
  canonical_sort_edge(ordered_edges_v1);

  for (unsigned int i = 0; i < edges_v1.size(); ++i) {
    if (!(edges_v1[i] == ordered_edges_v1[i])) {
      std::cout << "[" << i << "] different:\n";
      std::cout << edges_v1[i].get_smallest_vertex() << ", " << edges_v1[i].get_biggest_vertex() << ", "
                << edges_v1[i].get_filtration_value() << ", " << edges_v1[i].get_direction() << "\n";
      std::cout << ordered_edges_v1[i].get_smallest_vertex() << ", " << ordered_edges_v1[i].get_biggest_vertex() << ", "
                << ordered_edges_v1[i].get_filtration_value() << ", " << ordered_edges_v1[i].get_direction() << "\n";
    }
  }
}

void test_edges_asymetry(const std::vector<Point>& points, double nu, double mu, ORE::Order_policy p) {
  auto comp = [](const ZE& e1, const ZE& e2) {
    if (e1.get_smallest_vertex() != e2.get_smallest_vertex())
      return e1.get_smallest_vertex() < e2.get_smallest_vertex();
    return e1.get_biggest_vertex() < e2.get_biggest_vertex();
  };
  std::set<ZE, std::function<bool(ZE, ZE)> > curedges(comp);
  for (const auto& e : ORE::get_iterator_range(nu, mu, points, Gudhi::Euclidean_distance(), p)) {
    if (e.get_direction())
      curedges.insert(e);
    else
      curedges.erase(e);
  }
  std::cout << "final state: " << curedges.size() << "\n";
  for (const auto& e : curedges) {
    std::cout << e.get_smallest_vertex() << ", " << e.get_biggest_vertex() << ", " << e.get_filtration_value() << ", "
              << e.get_direction() << "\n";
  }
}

void test_edges_timings(const std::vector<Point>& points, double nu, double mu, ORE::Order_policy p) {
//   {
//     Gudhi::Clock time1("Vector version");
//     std::vector<ZE> edges_v1 = ORE::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);
//     std::cout << edges_v1.size() << "\n";
//     time1.end();
//     std::cout << time1;
//   }

  {
    Gudhi::Clock time2("Iterator version");
    unsigned int i = 0;
    for ([[maybe_unused]] const auto& e : ORE::get_iterator_range(nu, mu, points, Gudhi::Euclidean_distance(), p)) {
      ++i;
    }
    std::cout << i << "\n";
    time2.end();
    std::cout << time2;
  }

  {
    Gudhi::Clock time1("Vector version");
    std::vector<ZE> edges_v1 = ORE::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);
    std::cout << edges_v1.size() << "\n";
    time1.end();
    std::cout << time1;
  }
}

void test_edges(const std::vector<Point>& points, double nu, double mu) {
  // ORE::Order_policy p = ORE::Order_policy::FARTHEST_POINT_ORDERING;
  // ORE::Order_policy p = ORE::Order_policy::ALREADY_ORDERED;
  ORE::Order_policy p = ORE::Order_policy::RANDOM_POINT_ORDERING;

  // test_edges_comp(points, nu, mu, p);
  // test_edges_canonical_sort(points, nu, mu, p);
  // test_edges_asymetry(points, nu, mu, p);
  test_edges_timings(points, nu, mu, p);
}

void test_simplices_print(const std::vector<Point>& points, double nu, double mu, int maxDim, ORE::Order_policy p) {
  ST st;

  auto start = ORE::begin(nu, mu, points, Gudhi::Euclidean_distance(), p);
  auto end = ORE::end();
  for (auto& t : ORi::get_iterator_range(start, end, st, maxDim)) {
    for (auto v : st.simplex_vertex_range(std::get<0>(t))) std::cout << v << " ";
    std::cout << " -- " << std::get<1>(t) << ", " << std::get<2>(t) << "\n";
  }
}

void test_simplices_comp(const std::vector<Point>& points, double nu, double mu, int maxDim, ORE::Order_policy p){
  ST st;

  auto startEIt = ORE::begin(nu, mu, points, Gudhi::Euclidean_distance(), p);
  auto endEIt = ORE::end();
  auto vec = ORE::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);
  auto startEVec = vec.begin();
  auto endEVec = vec.end();

  auto startIt = ORi::begin(startEIt, endEIt, st, maxDim);
  auto endIt = ORi::end();
  auto rangeVec = ORv::get_iterator_range(startEVec, endEVec, st, maxDim);
  auto startVec = rangeVec.begin();
  auto endVec = rangeVec.end();
  for (; startIt != endIt && startVec != endVec; ++startIt,++startVec) {
    if ()
    for (auto v : st.simplex_vertex_range(std::get<0>(t))) std::cout << v << " ";
    std::cout << " -- " << std::get<1>(t) << ", " << std::get<2>(t) << "\n";
  }
}

void test_simplices(const std::vector<Point>& points, double nu, double mu, int maxDim) {
  ORE::Order_policy p = ORE::Order_policy::FARTHEST_POINT_ORDERING;
  // ORE::Order_policy p = ORE::Order_policy::ALREADY_ORDERED;
  // ORE::Order_policy p = ORE::Order_policy::RANDOM_POINT_ORDERING;
  
//   test_simplices_print(points, nu, mu, maxDim, p);
  test_simplices_comp(points, nu, mu, maxDim, p);
}

int main(int argc, char* const argv[]) {
  if (argc != 5 && argc != 6) {
    std::cout << "Usage: ./comp nu mu max_dim nomberOfPoints [seed]\n";
    return 0;
  }

  double nu = std::stod(argv[1]);
  double mu = std::stod(argv[2]);
  int maxDim = std::stoi(argv[3]);
  unsigned int numberOfPoints = std::stoi(argv[4]);
  int seed = -1;

  if (argc == 6) seed = std::stoi(argv[5]);

  std::cout << "nu, mu: " << nu << ", " << mu << "\n";
  std::cout << "number of points: " << numberOfPoints << "\n";
  std::cout << "seed: " << seed << "\n";

  std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

//   test_edges(points, nu, mu);
  test_simplices(points, nu, mu, maxDim);

  return 0;
}
