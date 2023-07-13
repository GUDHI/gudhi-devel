/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2014 Inria
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

#include <boost/range/adaptors.hpp>

#include <gudhi/Zigzag_persistence.h>
#include <gudhi/Simplex_tree.h>
#include <fzz/fzz.h>
#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/zigzag-persistence.h>
#include <dionysus/format.h>

#include "rips-zigzag-dionysus.h"

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_wide_indexation>;
using Gudhi::persistence_matrix::Zigzag_options;
using CT = Gudhi::persistence_matrix::Column_types;
using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<ST, Zigzag_options<CT::VECTOR> >;
using Vertex_handle = ST::Vertex_handle;
using Filtration_value = ST::Filtration_value;
using interval_filtration = ZP::interval_filtration;

using DField = dionysus::Z2Field;
using Simplex = dionysus::Simplex<>;
using DFiltration = dionysus::Filtration<Simplex>;
using DZZ = dionysus::ZigzagPersistence<DField>;
using DIndex = typename DZZ::Index;
using DChain = dionysus::ChainEntry<DField, typename DFiltration::Cell>;
using DIChain = dionysus::ChainEntry<DField, DIndex>;

void print_complex(ZP& zp){
	std::clog << std::endl << "Current complex:" << std::endl;
	zp.print_current_complex();
}

// void print_barcode(ZP& zp){
// 	std::clog << std::endl << "Current barcode:" << std::endl;
// 	for (auto& bar : zp.persistence_diagram()){
// 		std::clog << std::floor(bar.birth()) << " - ";
// 		if (bar.death() == std::numeric_limits<Filtration_value>::infinity()){
// 			std::clog << "inf";
// 		} else {
// 			std::clog << std::floor(bar.death());
// 		}
// 		std::clog << " (" << bar.dim() << ")\n";
// 	}
// }

std::vector< std::pair<unsigned int, unsigned int> > print_indices(ZP& zp, unsigned int numberOfSimplices){
	std::set<unsigned int> essentials;
	std::vector< std::pair<unsigned int, unsigned int> > res;

	for (unsigned int i = 0; i < numberOfSimplices; ++i){
		essentials.insert(essentials.end(), i);
	}

	for (auto& bar : zp.index_persistence_diagram()){
		// std::clog << bar.birth() << " - ";
		// std::clog << bar.death();
		// std::clog << " (" << bar.dim() << ")\n";
		res.emplace_back(bar.birth(), bar.death());
		essentials.erase(bar.birth());
		essentials.erase(bar.death());
	}

	for (unsigned int v : essentials){
		// std::clog << v << " - ";
		// std::clog << "inf\n";
		res.emplace_back(v, numberOfSimplices);
	}

	return res;
}

// std::vector<std::vector<Vertex_handle> > get_simplices()
// {
// 	return {
// 		{0}, 
// 		{1}, 
// 		{2}, 
// 		{0,1},
// 		{0,2},
// 		{3},
// 		{1,2},
// 		{4},
// 		{3,4},
// 		{5},
// 		{0,1,2},
// 		{4,5},
// 		{3,5},
// 		{3,4,5},
// 		{0,1,2},			//r
// 		{3,4,5},			//r
// 		{1,4},
// 		{0,1,2},
// 		{2,4},
// 		{3,4,5},
// 		{0,4},
// 		{0,2,4},
// 		{1,2,4},
// 		{0,1,4},
// 		{3,4,5},			//r
// 		{3,4},					//r
// 		{3,5},					//r
// 		{0,1,2,4}};
// }

// std::vector<Filtration_value> get_filtration_values()
// {
// 	return {
// 		0, 0, 0,
// 		1, 1, 1,
// 		2, 2, 2,
// 		3, 3, 3, 3,
// 		4,
// 		5,
// 		6, 6, 6,
// 		7, 7, 7, 7, 7, 7,
// 		8,
// 		9, 9, 9
// 	};
// }

// std::vector<bool> get_directions()
// {
// 	return {
// 		true, 
// 		true, 
// 		true, 
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		false,
// 		false,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		true,
// 		false,
// 		false,
// 		false,
// 		true};
// }

std::vector< std::pair<unsigned int, unsigned int> > compute_with_gudhi(
	const std::vector<std::vector<Vertex_handle> >& simplices,
	const std::vector<bool>& dirs)
{
	ZP zp(simplices.size());

	std::cout << "====================== Gudhi =====================\n";

	// std::vector<Filtration_value> filValues = get_filtration_values();
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
		// print_complex(zp);
	}

	std::cout << "==================================================\n";

	// print_complex(zp);
	// print_barcode(zp);
	return print_indices(zp, i);
}

std::vector< std::pair<unsigned int, unsigned int> > compute_with_dionysus(
	const std::vector<std::vector<Vertex_handle> >& simplices,
	const std::vector<bool>& dirs)
{
	DField k;
	std::unordered_map<Simplex, unsigned int> indices;
	DZZ persistence(k);
	std::vector< std::pair<unsigned int, unsigned int> > res;

	std::cout << "==================== Dionysus ====================\n";

	std::set<unsigned int> essentials;

	unsigned int op = 0;
	unsigned int idx = 0;

	for (const std::vector<Vertex_handle>& simplex : simplices){
		Simplex c(simplex);
		DIndex pair;
		if (dirs[op]) {
			indices.try_emplace(c, idx++);
			// int dim = boost::distance(c.boundary(persistence.field()));
			// dim = dim == 0 ? 0 : dim -1;
			// fmt::print("[{}] Adding: {} : {}\n", op, c, dim);
			pair = persistence.add(c.boundary(persistence.field()) |
								   boost::adaptors::transformed([&indices](const DChain& e) {
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

	for (unsigned int v : essentials){
		// fmt::print("{} - inf\n", v);
		res.emplace_back(v, op);
	}

	std::cout << "==================================================\n";

	return res;
}

std::vector< std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer> > compute_with_fzz(
	const std::vector<std::vector<Vertex_handle> >& simplices,
	const std::vector<bool>& dirs)
{
	std::vector< std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer> > persistence;
	// std::vector< std::pair<unsigned int, unsigned int> > res;
	FZZ::FastZigzag fzz;

	std::cout << "======================= FZZ ======================\n";

	fzz.compute(simplices, dirs, &persistence);

	std::sort(persistence.begin(), persistence.end(), 
			[](const std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer>& p1, const std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer>& p2){
				if (std::get<1>(p1) == std::get<1>(p2)){
					return std::get<0>(p1) < std::get<0>(p2);
				}

				return std::get<1>(p1) < std::get<1>(p2);
			});

	// for (auto& t : persistence)
	// 	res.emplace_back(std::get<0>(t), std::get<1>(t));

	// for (const auto& e : persistence) {
	// 	std::cout << (std::get<0>(e) - 1) << " - ";
	// 	if (static_cast<unsigned int>(std::get<1>(e)) == simplices.size()) std::cout << "inf";
	// 	else std::cout << std::get<1>(e);
	// 	std::cout << " (" << std::get<2>(e) << ")" << std::endl;
	// }

	std::cout << "==================================================\n";

	return persistence;
}

bool are_equal(const std::vector<std::pair<unsigned int, unsigned int> >& gudhiRes,
			   const std::vector<std::pair<unsigned int, unsigned int> >& dioRes) 
{
	if (gudhiRes.size() != dioRes.size()) return false;

	for (unsigned int i = 0; i < gudhiRes.size(); ++i){
		if (gudhiRes[i].first != dioRes[i].first || gudhiRes[i].second != dioRes[i].second)
			return false;
	}

	return true;
}

bool are_equal(const std::vector<std::pair<unsigned int, unsigned int> >& gudhiRes,
			   const std::vector<std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer> >& fzzRes) 
{
	if (gudhiRes.size() != fzzRes.size()) return false;

	for (unsigned int i = 0; i < gudhiRes.size(); ++i){
		if (static_cast<int>(gudhiRes[i].first) != std::get<0>(fzzRes[i]) - 1 || static_cast<int>(gudhiRes[i].second) != std::get<1>(fzzRes[i]))
			return false;
	}

	return true;
}

void print(const std::vector<std::pair<unsigned int, unsigned int> >& res, unsigned int infValue){
	for (const auto& p : res) {
		std::cout << p.first << " - ";
		if (p.second == infValue) std::cout << "inf";
		else std::cout << p.second;
		std::cout << std::endl;
	}
}

void print(const std::vector<std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer> >& res, int infValue){
	for (const auto& e : res) {
		std::cout << (std::get<0>(e) - 1) << " - ";
		if (std::get<1>(e) == infValue) std::cout << "inf";
		else std::cout << std::get<1>(e);
		std::cout << std::endl;
	}
}

void print_differences(const std::vector<std::pair<unsigned int, unsigned int> >& gudhiRes,
					   const std::vector<std::pair<unsigned int, unsigned int> >& dioRes, 
					   unsigned int infValue) 
{
	for (unsigned int i = 0; i < gudhiRes.size(); ++i){
		if (gudhiRes[i].first != dioRes[i].first || gudhiRes[i].second != dioRes[i].second){
			std::string dg = gudhiRes[i].second == infValue ? "inf" : std::to_string(gudhiRes[i].second);
			std::string dd = dioRes[i].second == infValue ? "inf" : std::to_string(dioRes[i].second);
			std::cout << "[" << i << "] " 
					  << gudhiRes[i].first << " - " << dg 
					  << " / " 
					  << dioRes[i].first << " - " << dd 
					  << "\n";
		}
	}
}

void print_differences(const std::vector<std::pair<unsigned int, unsigned int> >& gudhiRes,
					   const std::vector<std::tuple<FZZ::Integer, FZZ::Integer, FZZ::Integer> >& fzzRes, 
					   int infValue) 
{
	for (unsigned int i = 0; i < gudhiRes.size(); ++i){
		if (static_cast<int>(gudhiRes[i].first) != std::get<0>(fzzRes[i]) || static_cast<int>(gudhiRes[i].second) != std::get<1>(fzzRes[i])){
			std::string dg = static_cast<int>(gudhiRes[i].second) == infValue ? "inf" : std::to_string(gudhiRes[i].second);
			std::string dd = std::get<1>(fzzRes[i]) == infValue ? "inf" : std::to_string(std::get<1>(fzzRes[i]));
			std::cout << "[" << i << "] " 
					  << gudhiRes[i].first << " - " << dg 
					  << " / " 
					  << std::get<0>(fzzRes[i]) << " - " << dd 
					  << "\n";
		}
	}
}



int main(int argc, char* const argv[]) {
	if (argc < 2 || argc > 3) {
		std::cout << "Wrong arguments.\n";
		return 0;
	}

	// std::vector<std::vector<Vertex_handle> > simplices = get_simplices();
	// std::vector<bool> dirs = get_directions();
	std::vector<std::vector<Vertex_handle> > simplices;
	std::vector<bool> dirs;

	unsigned int numberOfPoints = std::stoi(argv[1]);
	int seed = -1;

	if (argc == 3)
		seed = std::stoi(argv[2]);

	unsigned int numberOfSimplices = build_rips_zigzag_filtration(simplices, dirs, numberOfPoints, seed);
	std::cout << "\n" << "numberOfSimplices: " << numberOfSimplices << "\n";
	
	auto gudhiRes = compute_with_gudhi(simplices, dirs);
	auto dioRes = compute_with_dionysus(simplices, dirs);
	auto fzzRes = compute_with_fzz(simplices, dirs);

	std::cout << "Res sizes: " << gudhiRes.size() << ", " << dioRes.size() << ", " << fzzRes.size() << "\n";

	bool firstRes = are_equal(gudhiRes, dioRes);
	if (!firstRes){
		std::cout << "------------------------ Gudhi and Dionysus results are not equal!\n";
		// print(gudhiRes, numberOfSimplices);
		// std::cout << "------------------------\n";
		// print(dioRes, numberOfSimplices);
		print_differences(gudhiRes, dioRes, numberOfSimplices);
	} else {
		std::cout << "+++++++++++++++++++++++++ Gudhi and Dionysus results are equal.\n";
	}

	if (!are_equal(gudhiRes, fzzRes)){
		std::cout << "------------------------ Gudhi and FZZ results are not equal!\n";
		// if (firstRes) {
		// 	print(gudhiRes, numberOfSimplices);
		// 	std::cout << "------------------------\n";
		// }
		// print(fzzRes, numberOfSimplices);
	} else {
		std::cout << "+++++++++++++++++++++++++ Gudhi and FZZ results are equal.\n";
	}

	return 0;
}
