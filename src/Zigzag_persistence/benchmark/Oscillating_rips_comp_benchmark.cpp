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
#include <iostream>
#include <string>
#include <utility>  // for pair
#include <vector>

#include <benchmark/benchmark.h>
#include <gudhi/Zigzag_persistence/oscillating_rips_iterators.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Clock.h>
#include <gudhi/Oscillating_rips_persistence.h>

#include "ext_zz/fzz/fzz.h"
#include "ext_zz/dionysus/simplex.h"
#include "ext_zz/dionysus/filtration.h"
#include "ext_zz/dionysus/zigzag-persistence.h"

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_oscillating_rips>;
using Vertex_handle = ST::Vertex_handle;
using Filtration_value = ST::Filtration_value;
using Simplex_handle = ST::Simplex_handle;
using Square = Gudhi::zigzag_persistence::Square_root_edge_modifier<Filtration_value>;
using ZE = Gudhi::zigzag_persistence::Zigzag_edge<Filtration_value>;
using ORE = Gudhi::zigzag_persistence::Oscillating_rips_edge_range<Filtration_value>;
using OR = Gudhi::zigzag_persistence::Oscillating_rips_simplex_range<ST,ORE::Oscillating_rips_edge_iterator>;
using ORv = Gudhi::zigzag_persistence::Oscillating_rips_simplex_range<ST,std::vector<ZE>::const_iterator>;
using Point = std::vector<double>;

using DField = dionysus::Z2Field;
using Simplex = dionysus::Simplex<>;
using DFiltration = dionysus::Filtration<Simplex>;
using DZZ = dionysus::ZigzagPersistence<DField>;
using DIndex = typename DZZ::Index;
using DChain = dionysus::ChainEntry<DField, typename DFiltration::Cell>;
using DIChain = dionysus::ChainEntry<DField, DIndex>;

struct Simplex_tree_options_test {
  typedef Gudhi::linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::int32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = true;
  static const bool link_nodes_by_label = false;
  static const bool simplex_handle_strong_validity = false;
};

// using ST_std = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_negative_indexation>;
using ST_std = Gudhi::Simplex_tree<Simplex_tree_options_test>;
using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<ST_std>;
using Barcode = std::vector<typename ZP::filtration_value_interval>;

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

  return finalPoints;
}

std::vector<std::pair<unsigned int, unsigned int> > compute_with_dionysus(const std::vector<ZE>& edges, int maxDim) {
  ST st;
  DField k;
  std::unordered_map<Simplex, unsigned int> indices;
  DZZ persistence(k);
  std::vector<std::pair<unsigned int, unsigned int> > res;

  std::set<unsigned int> essentials;

  unsigned int op = 0;
  unsigned int idx = 0;

  auto start = edges.begin();
  auto end = edges.end();
  for (const auto& t : ORv::get_iterator_range(start, end, st, maxDim)) {
    auto r = st.simplex_vertex_range(std::get<0>(t));
    std::vector<Vertex_handle> simplex(r.begin(), r.end());
    Simplex c(simplex);
    DIndex pair;
    if (std::get<2>(t)) {
      indices.try_emplace(c, idx++);
      // int dim = boost::distance(c.boundary(persistence.field()));
      // dim = dim == 0 ? 0 : dim -1;
      // fmt::print("[{}] Adding: {} : {}\n", op, c, dim);
      pair =
          persistence.add(c.boundary(persistence.field()) | boost::adaptors::transformed([&indices](const DChain& e) {
                            return DIChain(e.element(), indices.find(e.index())->second);
                          }));
    } else {
      // fmt::print("[{}] Removing: {} : {}\n", op, c, boost::distance(c.boundary(persistence.field())) - 1);
      auto idxIt = indices.find(c);
      pair = persistence.remove(idxIt->second);
      indices.erase(idxIt);
    }

    if (pair != DZZ::unpaired()) {
      // fmt::print("{} - {}\n", pair, op);
      res.emplace_back(pair, op);
      essentials.erase(pair);
    } else {
      essentials.insert(essentials.end(), op);
    }
    op++;
  }

  for (unsigned int v : essentials) {
    // fmt::print("{} - inf\n", v);
    res.emplace_back(v, op);
  }

  return res;
}

std::vector<std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer> > compute_with_fzz(const std::vector<ZE>& edges,
                                                                                    int maxDim) {
  std::vector<std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer> > persistence;
  FZZ::FastZigzag fzz;
  std::vector<std::vector<Vertex_handle> > simplices;
  std::vector<bool> dirs;
  ST st;

  auto start = edges.begin();
  auto end = edges.end();
  for (const auto& t : ORv::get_iterator_range(start, end, st, maxDim)){
    auto r = st.simplex_vertex_range(std::get<0>(t));
    simplices.emplace_back(r.begin(), r.end());
    dirs.push_back(std::get<2>(t));
  }

  fzz.compute(simplices, dirs, &persistence);

  std::sort(persistence.begin(), persistence.end(),
            [](const std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer>& p1,
               const std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer>& p2) {
              if (std::get<1>(p1) == std::get<1>(p2)) {
                return std::get<0>(p1) < std::get<0>(p2);
              }

              return std::get<1>(p1) < std::get<1>(p2);
            });

  return persistence;
}

Barcode compute_with_gudhi_v2(const std::vector<Point>& points, double nu, double mu, int maxDim) {
  ZP zp;
  ST st;

  auto start = ORE::begin(nu, mu, points, Gudhi::Euclidean_distance(), ORE::Order_policy::FARTHEST_POINT_ORDERING);
  auto end = ORE::end();
  for (const auto& t : OR::get_iterator_range(start, end, st, maxDim)) {
    auto r = st.simplex_vertex_range(std::get<0>(t));
    std::vector<Vertex_handle> simplex(r.begin(), r.end());
    if (std::get<2>(t))
      zp.insert_simplex(simplex, std::get<1>(t));
    else
      zp.remove_simplex(simplex, std::get<1>(t));
  }

  return zp.get_persistence_diagram();
}

static void ORP_with_gudhi(benchmark::State& state) {
  double nu = state.range(0);
  double mu = state.range(1);
  int maxDim = state.range(2);
  unsigned int numberOfPoints = state.range(3);
  int seed = 0;

  std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

  for (auto _ : state) {
    auto res = Gudhi::zigzag_persistence::compute_oscillating_rips_persistence(points, nu, mu, maxDim);
  }
}
BENCHMARK(ORP_with_gudhi)
    ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
                   benchmark::CreateDenseRange(3, 5, 1),
                   benchmark::CreateRange(2, 5, 4),
                   benchmark::CreateRange(5, 5, 4)})
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(ORP_with_gudhi)
    ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
                   benchmark::CreateDenseRange(3, 5, 1),
                   benchmark::CreateRange(4, 50, 4),
                   benchmark::CreateRange(20, 64, 4)})
    ->Unit(benchmark::kMillisecond);
// BENCHMARK(ORP_with_gudhi)
//     ->ArgsProduct({benchmark::CreateDenseRange(1, 2, 1), 
//                    benchmark::CreateDenseRange(3, 4, 1),
//                    benchmark::CreateRange(4, 50, 4),
//                    benchmark::CreateRange(256, 500, 4)})
//     ->Unit(benchmark::kMillisecond);

// static void ORP_with_gudhi_v2(benchmark::State& state) {
//   double nu = state.range(0);
//   double mu = state.range(1);
//   int maxDim = state.range(2);
//   unsigned int numberOfPoints = state.range(3);
//   int seed = 0;

//   std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

//   for (auto _ : state) {
//     auto res = compute_with_gudhi_v2(points, nu, mu, maxDim);
//   }
// }
// BENCHMARK(ORP_with_gudhi_v2)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
//                    benchmark::CreateDenseRange(3, 5, 1),
//                    benchmark::CreateRange(2, 5, 4),
//                    benchmark::CreateRange(5, 5, 4)})
//     ->Unit(benchmark::kMicrosecond);
// BENCHMARK(ORP_with_gudhi_v2)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
//                    benchmark::CreateDenseRange(3, 5, 1),
//                    benchmark::CreateRange(4, 50, 4),
//                    benchmark::CreateRange(20, 64, 4)})
//     ->Unit(benchmark::kMillisecond);
// // BENCHMARK(ORP_with_gudhi_v2)
// //     ->ArgsProduct({benchmark::CreateDenseRange(1, 2, 1), 
// //                    benchmark::CreateDenseRange(3, 4, 1),
// //                    benchmark::CreateRange(4, 50, 4),
// //                    benchmark::CreateRange(256, 500, 4)})
// //     ->Unit(benchmark::kMillisecond);

static void ORP_with_dionysus(benchmark::State& state) {
  double nu = state.range(0);
  double mu = state.range(1);
  int maxDim = state.range(2);
  unsigned int numberOfPoints = state.range(3);
  int seed = 0;

  std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

  for (auto _ : state) {
    auto vec = ORE::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance());
    auto res = compute_with_dionysus(vec, maxDim);
  }
}
BENCHMARK(ORP_with_dionysus)
    ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
                   benchmark::CreateDenseRange(3, 5, 1),
                   benchmark::CreateRange(2, 5, 4),
                   benchmark::CreateRange(5, 5, 4)})
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(ORP_with_dionysus)
    ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
                   benchmark::CreateDenseRange(3, 5, 1),
                   benchmark::CreateRange(4, 50, 4),
                   benchmark::CreateRange(20, 64, 4)})
    ->Unit(benchmark::kMillisecond);
// BENCHMARK(ORP_with_dionysus)
//     ->ArgsProduct({benchmark::CreateDenseRange(1, 2, 1), 
//                    benchmark::CreateDenseRange(3, 4, 1),
//                    benchmark::CreateRange(4, 50, 4),
//                    benchmark::CreateRange(256, 500, 4)})
//     ->Unit(benchmark::kMillisecond);

static void ORP_with_fzz(benchmark::State& state) {
  double nu = state.range(0);
  double mu = state.range(1);
  int maxDim = state.range(2);
  unsigned int numberOfPoints = state.range(3);
  int seed = 0;

  std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

  for (auto _ : state) {
    auto vec = ORE::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance());
    auto res = compute_with_fzz(vec, maxDim);
  }
}
BENCHMARK(ORP_with_fzz)
    ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
                   benchmark::CreateDenseRange(3, 5, 1),
                   benchmark::CreateRange(2, 5, 4),
                   benchmark::CreateRange(5, 5, 4)})
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(ORP_with_fzz)
    ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
                   benchmark::CreateDenseRange(3, 5, 1),
                   benchmark::CreateRange(4, 50, 4),
                   benchmark::CreateRange(20, 64, 4)})
    ->Unit(benchmark::kMillisecond);
// BENCHMARK(ORP_with_fzz)
//     ->ArgsProduct({benchmark::CreateDenseRange(1, 2, 1), 
//                    benchmark::CreateDenseRange(3, 4, 1),
//                    benchmark::CreateRange(4, 50, 4),
//                    benchmark::CreateRange(256, 500, 4)})
//     ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();

