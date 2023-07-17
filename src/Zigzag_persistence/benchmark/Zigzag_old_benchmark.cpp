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

#include <gudhi/Zigzag_persistence_old.h>
#include <gudhi/Simplex_tree.h>

#include "rips-zigzag-dionysus.h"

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_wide_indexation>;
using coltype = Gudhi::zigzag_persistence::Zigzag_persistence_collist;
using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<ST, coltype>;
// using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<ST>;
using Vertex_handle = ST::Vertex_handle;
using Filtration_value = ST::Filtration_value;
using interval_filtration = ZP::interval_filtration;

std::vector< std::pair<unsigned int, unsigned int> > print_indices(ZP& zp, unsigned int numberOfSimplices){
	std::set<unsigned int> essentials;
	std::vector< std::pair<unsigned int, unsigned int> > res;

	for (unsigned int i = 0; i < numberOfSimplices; ++i){
		essentials.insert(essentials.end(), i);
	}

	for (auto& bar : zp.index_persistence_diagram()){
		res.emplace_back(bar.birth(), bar.death());
		essentials.erase(bar.birth());
		essentials.erase(bar.death());
	}

	for (unsigned int v : essentials){
		res.emplace_back(v, numberOfSimplices);
	}

	return res;
}

std::vector< std::pair<unsigned int, unsigned int> > compute_with_gudhi(
	const std::vector<std::vector<Vertex_handle> >& simplices,
	const std::vector<bool>& dirs)
{
	ZP zp(simplices.size());

	// std::cout << "====================== Gudhi =====================\n";

	std::vector<Filtration_value> filValues(simplices.size(), 1.0);

	auto start = simplices.begin();
	auto filIt = filValues.begin();
	unsigned int i = 0;

	while (start != simplices.end()){
		unsigned int c = 1;
		auto end = start + 1;
		++i;
		while (end != simplices.end() && dirs[i - 1] == dirs[i]) {
			++end;
			++i;
			++c;
		}

		if (dirs[i - 1]){
			zp.insert_simplices_contiguously(
				start, end, filIt);
		} else {
			zp.remove_simplices_contiguously(
				start, end, filIt);
		}

		start = end;
		filIt += c;
	}

	// std::cout << "==================================================\n";

	return print_indices(zp, i);
}

static void Compute_zigzag_with_gudhi(benchmark::State& state) {
	unsigned int numberOfPoints = state.range(0);;
	int seed = numberOfPoints;
	std::vector<std::vector<Vertex_handle> > simplices;
	std::vector<bool> dirs;

	build_rips_zigzag_filtration(simplices, dirs, numberOfPoints, seed);

	for (auto _ : state){
		compute_with_gudhi(simplices, dirs);
	}
}
BENCHMARK(Compute_zigzag_with_gudhi)->RangeMultiplier(2)->Range(100, 1000)->Unit(benchmark::kMillisecond)->MinWarmUpTime(1);

BENCHMARK_MAIN();

