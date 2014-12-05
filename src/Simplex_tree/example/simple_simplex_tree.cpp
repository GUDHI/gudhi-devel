/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <ctime>
#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

typedef std::vector< Vertex_handle > typeVectorVertex;
typedef std::pair<typeVectorVertex, Filtration_value> typeSimplex;
typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;

int main (int argc, char * const argv[])
{
	const Filtration_value FIRST_FILTRATION_VALUE = 0.1;
	const Filtration_value SECOND_FILTRATION_VALUE = 0.2;
	const Filtration_value THIRD_FILTRATION_VALUE = 0.3;
	const Filtration_value FOURTH_FILTRATION_VALUE = 0.4;
	Vertex_handle FIRST_VERTEX_HANDLE =  (Vertex_handle)0;
	Vertex_handle SECOND_VERTEX_HANDLE = (Vertex_handle)1;
	Vertex_handle THIRD_VERTEX_HANDLE = (Vertex_handle)2;
	Vertex_handle FOURTH_VERTEX_HANDLE = (Vertex_handle)3;

	// TEST OF INSERTION
	std::cout << "********************************************************************" << std::endl;
	std::cout << "EXAMPLE OF SIMPLE INSERTION" << std::endl;
	//Construct the Simplex Tree
	Simplex_tree<> simplexTree;

	// Simplex to be inserted:
	//    1
	//    o
	//   /X\
	//  o---o---o
	//  2   0   3

	// ++ FIRST
	std::cout << "   * INSERT 0" << std::endl;
	typeVectorVertex firstSimplexVector;
	firstSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	typeSimplex firstSimplex = std::make_pair(firstSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
	typePairSimplexBool returnValue =
			simplexTree.insert ( firstSimplex.first, firstSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + 0 INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - 0 NOT INSERTED" << std::endl;
	}

	// ++ SECOND
	std::cout << "   * INSERT 1" << std::endl;
	typeVectorVertex secondSimplexVector;
	secondSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	typeSimplex secondSimplex = std::make_pair(secondSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
	returnValue =
			simplexTree.insert ( secondSimplex.first, secondSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + 1 INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - 1 NOT INSERTED" << std::endl;
	}

	// ++ THIRD
	std::cout << "   * INSERT (0,1)" << std::endl;
	typeVectorVertex thirdSimplexVector;
	thirdSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	thirdSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	typeSimplex thirdSimplex = std::make_pair(thirdSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
	returnValue =
			simplexTree.insert ( thirdSimplex.first, thirdSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + (0,1) INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - (0,1) NOT INSERTED" << std::endl;
	}

	// ++ FOURTH
	std::cout << "   * INSERT 2" << std::endl;
	typeVectorVertex fourthSimplexVector;
	fourthSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	typeSimplex fourthSimplex = std::make_pair(fourthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
	returnValue =
			simplexTree.insert ( fourthSimplex.first, fourthSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + 2 INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - 2 NOT INSERTED" << std::endl;
	}

	// ++ FIFTH
	std::cout << "   * INSERT (2,0)" << std::endl;
	typeVectorVertex fifthSimplexVector;
	fifthSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	fifthSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	typeSimplex fifthSimplex = std::make_pair(fifthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
	returnValue =
			simplexTree.insert ( fifthSimplex.first, fifthSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + (2,0) INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - (2,0) NOT INSERTED" << std::endl;
	}

	// ++ SIXTH
	std::cout << "   * INSERT (2,1)" << std::endl;
	typeVectorVertex sixthSimplexVector;
	sixthSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	sixthSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	typeSimplex sixthSimplex = std::make_pair(sixthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
	returnValue =
			simplexTree.insert ( sixthSimplex.first, sixthSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + (2,1) INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - (2,1) NOT INSERTED" << std::endl;
	}

	// ++ SEVENTH
	std::cout << "   * INSERT (2,1,0)" << std::endl;
	typeVectorVertex seventhSimplexVector;
	seventhSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	seventhSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	seventhSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	typeSimplex seventhSimplex = std::make_pair(seventhSimplexVector, Filtration_value(THIRD_FILTRATION_VALUE));
	returnValue =
			simplexTree.insert ( seventhSimplex.first, seventhSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + (2,1,0) INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - (2,1,0) NOT INSERTED" << std::endl;
	}

	// ++ EIGHTH
	std::cout << "   * INSERT 3" << std::endl;
	typeVectorVertex eighthSimplexVector;
	eighthSimplexVector.push_back(FOURTH_VERTEX_HANDLE);
	typeSimplex eighthSimplex = std::make_pair(eighthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
	returnValue =
			simplexTree.insert ( eighthSimplex.first, eighthSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + 3 INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - 3 NOT INSERTED" << std::endl;
	}

	// ++ NINETH
	std::cout << "   * INSERT (3,0)" << std::endl;
	typeVectorVertex ninethSimplexVector;
	ninethSimplexVector.push_back(FOURTH_VERTEX_HANDLE);
	ninethSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	typeSimplex ninethSimplex = std::make_pair(ninethSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
	returnValue =
			simplexTree.insert ( ninethSimplex.first, ninethSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + (3,0) INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - (3,0) NOT INSERTED" << std::endl;
	}

	// ++ TENTH
	std::cout << "   * INSERT 0 (already inserted)" << std::endl;
	typeVectorVertex tenthSimplexVector;
	tenthSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	typeSimplex tenthSimplex = std::make_pair(tenthSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE)); // With a different filtration value
	returnValue =
			simplexTree.insert ( tenthSimplex.first, tenthSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + 0 INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - 0 NOT INSERTED" << std::endl;
	}

	// ++ ELEVENTH
	std::cout << "   * INSERT (2,1,0) (already inserted)" << std::endl;
	typeVectorVertex eleventhSimplexVector;
	eleventhSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	eleventhSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	eleventhSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	typeSimplex eleventhSimplex = std::make_pair(eleventhSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE));
	returnValue =
			simplexTree.insert ( eleventhSimplex.first, eleventhSimplex.second );

	if (returnValue.second == true)
	{
		std::cout << "   + (2,1,0) INSERTED" << std::endl;
		int nb_simplices = simplexTree.num_simplices() + 1;
		simplexTree.set_num_simplices(nb_simplices);
	}
	else
	{
		std::cout << "   - (2,1,0) NOT INSERTED" << std::endl;
	}

	// ++ GENERAL VARIABLE SET
	simplexTree.set_filtration(FOURTH_FILTRATION_VALUE); // Max filtration value
	simplexTree.set_dimension(2); // Max dimension = 2 -> (2,1,0)

	std::cout << "********************************************************************" << std::endl;
	// Display the Simplex_tree - Can not be done in the middle of 2 inserts
	std::cout << "* The complex contains " << simplexTree.num_simplices() << " simplices" << std::endl;
	std::cout << "   - dimension " << simplexTree.dimension() << "   - filtration " << simplexTree.filtration() << std::endl;
	std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
	for( auto f_simplex : simplexTree.filtration_simplex_range() )
	{
		std::cout << "   " << "[" << simplexTree.filtration(f_simplex) << "] ";
		for( auto vertex : simplexTree.simplex_vertex_range(f_simplex) )
		{
			std::cout << (int)vertex << " ";
		}
		std::cout << std::endl;
	}
	//   [0.1] 0
	//   [0.1] 1
	//   [0.1] 2
	//   [0.1] 3
	//   [0.2] 1 0
	//   [0.2] 2 0
	//   [0.2] 2 1
	//   [0.2] 3 0
	//   [0.3] 2 1 0
}
