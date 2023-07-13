/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <stdlib.h>    // atoi, getenv
#include <stddef.h>    // size_t
#include <stdio.h>     // printf
#include <string.h>    // strcmp
#include <unistd.h>    // read, write
#include <sys/stat.h>
#include <fcntl.h>

#include <gudhi/Zigzag_persistence.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Clock.h>

#include <iostream>
#include <utility>  // for pair
#include <vector>

#include "rips-zigzag-dionysus.h"

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_wide_indexation>;
using Gudhi::persistence_matrix::Zigzag_options;
using CT = Gudhi::persistence_matrix::Column_types;
using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<ST, Zigzag_options<CT::INTRUSIVE_LIST> >;
// using coltype = Gudhi::zigzag_persistence::Zigzag_persistence_collist;
// using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<ST, coltype>;
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

int main(int argc, char* const argv[]) {
	if (argc < 2 || argc > 3) {
		std::cout << "Wrong arguments.\n";
		return 0;
	}

	int perf_ctl_fd = open("/tmp/perf_ctl.fifo",O_WRONLY);
	int perf_ctl_ack_fd = open("/tmp/perf_ctl_ack.fifo",O_RDONLY);
	char ack[5];
	std::cout << "perf_ctl_fd: " << perf_ctl_fd << "\n";
	std::cout << "perf_ctl_ack_fd: " << perf_ctl_ack_fd << "\n";

	unsigned int numberOfPoints = std::stoi(argv[1]);
	int seed = -1;
	if (argc == 3)
		seed = std::stoi(argv[2]);
	
	std::vector<std::vector<Vertex_handle> > simplices;
	std::vector<bool> dirs;

	unsigned int numberOfSimplices = build_rips_zigzag_filtration(simplices, dirs, numberOfPoints, seed);
	std::cout << "Filtration size: " << simplices.size() << "\n";
	std::cout << "Number of simplices: " << numberOfSimplices << "\n";
	
	// Start the performance counter and read the ack
	if (perf_ctl_fd != -1){
		write(perf_ctl_fd, "enable\n", 8);
		read(perf_ctl_ack_fd, ack, 5);
		if(std::strcmp(ack, "ack\n") != 0){
			std::cout << "No acknowledgment\n";
			return 1;
		}
	}

	Gudhi::Clock time("Zigzag Rips");
	/* auto res =  */compute_with_gudhi(simplices, dirs);
	time.end();
	std::cout << time;

	// Stop the performance counter and read the ack
	if (perf_ctl_fd != -1){
		write(perf_ctl_fd, "disable\n", 9);
		read(perf_ctl_ack_fd, ack, 5);
		if(std::strcmp(ack, "ack\n") != 0){
			std::cout << "No acknowledgment\n";
			return 1;
		}
	}

	// for (const auto& p : res) {
	// 	std::cout << p.first << " - ";
	// 	if (p.second == numberOfSimplices) std::cout << "inf";
	// 	else std::cout << p.second;
	// 	std::cout << std::endl;
	// }

	return 0;
}
