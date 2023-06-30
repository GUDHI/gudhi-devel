/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>  // for pair
#include <vector>
#include <unordered_map>

#include <benchmark/benchmark.h>
#include <boost/range/adaptors.hpp>

#include <fzz/fzz.h>
#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/zigzag-persistence.h>
#include <dionysus/format.h>

#include "rips-zigzag-dionysus.h"

using DField = dionysus::Z2Field;
using Simplex = dionysus::Simplex<>;
using DFiltration = dionysus::Filtration<Simplex>;
using DZZ = dionysus::ZigzagPersistence<DField>;
using DIndex = typename DZZ::Index;
using DChain = dionysus::ChainEntry<DField, typename DFiltration::Cell>;
using DIChain = dionysus::ChainEntry<DField, DIndex>;
using Vertex_handle = int;

std::vector< std::pair<unsigned int, unsigned int> > compute_with_dionysus(
	const std::vector<std::vector<Vertex_handle> >& simplices,
	const std::vector<bool>& dirs)
{
	DField k;
	std::unordered_map<Simplex, unsigned int> indices;
	DZZ persistence(k);
	std::vector< std::pair<unsigned int, unsigned int> > res;

	// std::cout << "==================== Dionysus ====================\n";

	std::set<unsigned int> essentials;

	unsigned int op = 0;
	unsigned int idx = 0;

	for (const std::vector<Vertex_handle>& simplex : simplices){
		Simplex c(simplex);
		DIndex pair;
		if (dirs[op]) {
			indices.emplace(c, idx++);
			pair = persistence.add(c.boundary(persistence.field()) |
								   boost::adaptors::transformed([&indices](const DChain& e) {
										return DIChain(e.element(), indices.find(e.index())->second);
								   }));
		} else {
			auto idxIt = indices.find(c);
			pair = persistence.remove(idxIt->second);
			indices.erase(idxIt);
		}

		if (pair != DZZ::unpaired()) {
			res.emplace_back(pair, op);
			essentials.erase(pair);
		} else {
			essentials.insert(essentials.end(), op);
		}
		op++;
	}

	for (unsigned int v : essentials){
		res.emplace_back(v, op);
	}

	// std::cout << "==================================================\n";

	return res;
}

std::vector< std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer> > compute_with_fzz(
	const std::vector<std::vector<Vertex_handle> >& simplices,
	const std::vector<bool>& dirs)
{
	std::vector< std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer> > persistence;
	FZZ::FastZigzag fzz;

	// std::cout << "======================= FZZ ======================\n";

	fzz.compute(simplices, dirs, &persistence);

	std::sort(persistence.begin(), persistence.end(), 
			[](const std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer>& p1, const std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer>& p2){
				if (std::get<1>(p1) == std::get<1>(p2)){
					return std::get<0>(p1) < std::get<0>(p2);
				}

				return std::get<1>(p1) < std::get<1>(p2);
			});

	// std::cout << "==================================================\n";

	return persistence;
}

static void Compute_zigzag_with_dionysus(benchmark::State& state) {
	unsigned int numberOfPoints = state.range(0);;
	int seed = numberOfPoints;
	std::vector<std::vector<Vertex_handle> > simplices;
	std::vector<bool> dirs;

	build_rips_zigzag_filtration(simplices, dirs, numberOfPoints, seed);

	for (auto _ : state){
		compute_with_dionysus(simplices, dirs);
	}
}
BENCHMARK(Compute_zigzag_with_dionysus)->RangeMultiplier(2)->Range(100, 1000)->Unit(benchmark::kMillisecond)->MinWarmUpTime(1)->MinTime(1);

static void Compute_zigzag_with_fzz(benchmark::State& state) {
	unsigned int numberOfPoints = state.range(0);;
	int seed = numberOfPoints;
	std::vector<std::vector<Vertex_handle> > simplices;
	std::vector<bool> dirs;

	build_rips_zigzag_filtration(simplices, dirs, numberOfPoints, seed);

	for (auto _ : state){
		compute_with_fzz(simplices, dirs);
	}
}
BENCHMARK(Compute_zigzag_with_fzz)->RangeMultiplier(2)->Range(100, 1000)->Unit(benchmark::kMillisecond)->MinWarmUpTime(1)->MinTime(1);

BENCHMARK_MAIN();

