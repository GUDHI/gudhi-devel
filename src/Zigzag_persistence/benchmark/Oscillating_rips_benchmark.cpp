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

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_oscillating_rips>;
using Filtration_value = ST::Filtration_value;
using Simplex_handle = ST::Simplex_handle;
using Square = Gudhi::zigzag_persistence::Square_root_edge_modifier<Filtration_value>;
using ZE = Gudhi::zigzag_persistence::Zigzag_edge<Filtration_value>;
using ORE = Gudhi::zigzag_persistence::Oscillating_rips_edge_range<Filtration_value>;
using OR = Gudhi::zigzag_persistence::Oscillating_rips_simplex_range<ST,ORE::Oscillating_rips_edge_iterator>;
using ORv = Gudhi::zigzag_persistence::Oscillating_rips_simplex_range<ST,std::vector<ZE>::iterator>;
using Point = std::vector<double>;

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

// static void ORP_with_custom_iterator(benchmark::State& state) {
//   double nu = state.range(0);
//   double mu = state.range(1);
//   int maxDim = state.range(2);
//   unsigned int numberOfPoints = state.range(3);
//   int seed = 0;

//   std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

//   for (auto _ : state) {
//     auto res = Gudhi::zigzag_persistence::compute_oscillating_rips_persistence(points, nu, mu, maxDim, ORE::Order_policy::FARTHEST_POINT_ORDERING);
//   }
// }
// BENCHMARK(ORP_with_custom_iterator)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
//                    benchmark::CreateDenseRange(3, 5, 1),
//                    benchmark::CreateRange(2, 5, 4),
//                    benchmark::CreateRange(5, 5, 4)})
//     ->Unit(benchmark::kMicrosecond);
// BENCHMARK(ORP_with_custom_iterator)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
//                    benchmark::CreateDenseRange(3, 5, 1),
//                    benchmark::CreateRange(4, 50, 4),
//                    benchmark::CreateRange(20, 64, 4)})
//     ->Unit(benchmark::kMillisecond);
// BENCHMARK(ORP_with_custom_iterator)
//     ->ArgsProduct({benchmark::CreateDenseRange(1, 2, 1), 
//                    benchmark::CreateDenseRange(3, 4, 1),
//                    benchmark::CreateRange(4, 50, 4),
//                    benchmark::CreateRange(256, 500, 4)})
//     ->Unit(benchmark::kMillisecond);

// static void ORP_with_vector_iterator(benchmark::State& state) {
//   double nu = state.range(0);
//   double mu = state.range(1);
//   int maxDim = state.range(2);
//   unsigned int numberOfPoints = state.range(3);
//   int seed = 0;

//   std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

//   for (auto _ : state) {
//     auto res = Gudhi::zigzag_persistence::compute_oscillating_rips_persistence<
//         std::vector<Point>, Gudhi::zigzag_persistence::Edge_range_type::VECTOR>(points, nu, mu, maxDim, ORE::Order_policy::FARTHEST_POINT_ORDERING);
//   }
// }
// BENCHMARK(ORP_with_vector_iterator)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1),
//                    benchmark::CreateDenseRange(3, 5, 1),
//                    benchmark::CreateRange(2, 5, 4),
//                    benchmark::CreateRange(5, 5, 4)})
//     ->Unit(benchmark::kMicrosecond);
// BENCHMARK(ORP_with_vector_iterator)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1),
//                    benchmark::CreateDenseRange(3, 5, 1),
//                    benchmark::CreateRange(4, 50, 4),
//                    benchmark::CreateRange(20, 64, 4)})
//     ->Unit(benchmark::kMillisecond);
// BENCHMARK(ORP_with_vector_iterator)
//     ->ArgsProduct({benchmark::CreateDenseRange(1, 2, 1),
//                    benchmark::CreateDenseRange(3, 4, 1),
//                    benchmark::CreateRange(4, 50, 4),
//                    benchmark::CreateRange(256, 500, 4)})
//     ->Unit(benchmark::kMillisecond);

// static void ORP_with_custom_iterator_random(benchmark::State& state) {
//   double nu = state.range(0);
//   double mu = state.range(1);
//   int maxDim = state.range(2);
//   unsigned int numberOfPoints = state.range(3);
//   int seed = 0;

//   std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

//   for (auto _ : state) {
//     auto res = Gudhi::zigzag_persistence::compute_oscillating_rips_persistence(points, nu, mu, maxDim, ORE::Order_policy::RANDOM_POINT_ORDERING);
//   }
// }
// BENCHMARK(ORP_with_custom_iterator_random)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
//                    benchmark::CreateDenseRange(3, 5, 1),
//                    benchmark::CreateRange(4, 50, 4),     
//                    benchmark::CreateDenseRange(5, 15, 5)})
//     ->Unit(benchmark::kMillisecond);
// BENCHMARK(ORP_with_custom_iterator_random)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
//                    benchmark::CreateDenseRange(3, 4, 1),
//                    benchmark::CreateRange(4, 50, 4),     
//                    benchmark::CreateDenseRange(20, 50, 30)})
//     ->Unit(benchmark::kMillisecond);

static void ORP_with_vector_iterator_random(benchmark::State& state) {
  double nu = state.range(0);
  double mu = state.range(1);
  int maxDim = state.range(2);
  unsigned int numberOfPoints = state.range(3);
  int seed = 0;

  std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

  for (auto _ : state) {
    auto res = Gudhi::zigzag_persistence::compute_oscillating_rips_persistence<
        std::vector<Point>, Gudhi::zigzag_persistence::Edge_range_type::VECTOR>(points, nu, mu, maxDim, ORE::Order_policy::RANDOM_POINT_ORDERING);
  }
}
BENCHMARK(ORP_with_vector_iterator_random)
    ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
                   benchmark::CreateDenseRange(3, 5, 1),
                   benchmark::CreateRange(4, 50, 4),     
                   benchmark::CreateDenseRange(5, 15, 5)})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(ORP_with_vector_iterator_random)
    ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
                   benchmark::CreateDenseRange(3, 4, 1),
                   benchmark::CreateRange(4, 50, 4),     
                   benchmark::CreateDenseRange(20, 50, 30)})
    ->Unit(benchmark::kMillisecond);

// static void ORP_with_custom_iterator_ordered(benchmark::State& state) {
//   double nu = state.range(0);
//   double mu = state.range(1);
//   int maxDim = state.range(2);
//   unsigned int numberOfPoints = state.range(3);
//   int seed = 0;

//   std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

//   for (auto _ : state) {
//     auto res = Gudhi::zigzag_persistence::compute_oscillating_rips_persistence(points, nu, mu, maxDim, ORE::Order_policy::ALREADY_ORDERED);
//   }
// }
// BENCHMARK(ORP_with_custom_iterator_ordered)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
//                    benchmark::CreateDenseRange(3, 5, 1),
//                    benchmark::CreateRange(4, 50, 4),     
//                    benchmark::CreateRange(5, 50, 4)})
//     ->Unit(benchmark::kMillisecond);

// static void ORP_with_vector_iterator_ordered(benchmark::State& state) {
//   double nu = state.range(0);
//   double mu = state.range(1);
//   int maxDim = state.range(2);
//   unsigned int numberOfPoints = state.range(3);
//   int seed = 0;

//   std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

//   for (auto _ : state) {
//     auto res = Gudhi::zigzag_persistence::compute_oscillating_rips_persistence<
//         std::vector<Point>, Gudhi::zigzag_persistence::Edge_range_type::VECTOR>(points, nu, mu, maxDim, ORE::Order_policy::ALREADY_ORDERED);
//   }
// }
// BENCHMARK(ORP_with_vector_iterator_ordered)
//     ->ArgsProduct({benchmark::CreateDenseRange(0, 2, 1), 
//                    benchmark::CreateDenseRange(3, 5, 1),
//                    benchmark::CreateRange(4, 50, 4),     
//                    benchmark::CreateRange(5, 50, 4)})
//     ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();

