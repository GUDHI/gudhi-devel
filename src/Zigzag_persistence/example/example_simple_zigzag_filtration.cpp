/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Zigzag_persistence.h>
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <utility>  // for pair
#include <vector>

using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_wide_indexation>;
using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<Simplex_tree>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using interval_filtration = ZP::interval_filtration;

void print_complex(ZP& zp){
	std::clog << std::endl << "Current complex:" << std::endl;
	zp.print_current_complex();
}

void print_barcode(ZP& zp){
	std::clog << std::endl << "Current barcode:" << std::endl;
	for (auto& bar : zp.persistence_diagram()){
		std::clog << std::floor(bar.birth()) << " - ";
		if (bar.death() == std::numeric_limits<Filtration_value>::infinity()){
			std::clog << "inf";
		} else {
			std::clog << std::floor(bar.death());
		}
		std::clog << " (" << bar.dim() << ")\n";
	}
}

int main(int argc, char* const argv[]) {
	ZP zp;

	std::vector<std::vector<Vertex_handle> > simplices{
		{0},{1},{2},
		{0,1},{0,2},{3},
		{1,2},{4},{3,4},
		{5},{0,1,2},{4,5},{3,5}};
	std::vector<Filtration_value> fils{0,0,0,1,1,1,2,2,2,3,3,3,3};
	zp.insert_simplices_contiguously(simplices, fils);

	print_complex(zp);
	print_barcode(zp);

	std::vector<Vertex_handle> simplex{3,4,5};
	zp.insert_simplex(simplex, 4);

	print_complex(zp);
	print_barcode(zp);

	simplex[0] = 0;
	simplex[1] = 1;
	simplex[2] = 2;
	zp.remove_simplex(simplex, 5);

	print_complex(zp);
	print_barcode(zp);

	simplex[0] = 3;
	simplex[1] = 4;
	simplex[2] = 5;
	zp.remove_simplex(simplex, 6);

	print_complex(zp);
	print_barcode(zp);

	simplices = {{1,4},{0,1,2},{2,4},{3,4,5},{0,4},{0,2,4},{1,2,4},{0,1,4}};
	fils = {6,6,7,7,7,7,7,7};
	zp.insert_simplices_contiguously(simplices, fils);

	print_complex(zp);
	print_barcode(zp);

	simplices = {{3,4,5},{3,4},{3,5}};
	fils = {8,8,8};
	zp.remove_simplices_contiguously(simplices, fils);

	print_complex(zp);
	print_barcode(zp);

	simplex[0] = 0;
	simplex[1] = 1;
	simplex[2] = 2;
	simplex.push_back(4);
	zp.insert_simplex(simplex, 8);

	print_complex(zp);
	print_barcode(zp);

	return 0;
}
