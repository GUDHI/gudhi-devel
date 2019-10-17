/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Divyansh Pareek
 *
 *    Copyright (C) 2017 INRIA Sophia Antipolis (France)
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

#include <gudhi/MsMatrix.h>

int main()
{
	Simplex_tree stree;

	// Forming toy Simplicial Complex for example

	/* Complex to build. */

	/*    o 0            */
	/*    |              */
	/*    o 1            */
	/*   /"\             */
	/*  o---o            */
	/* 2 \X/ 3           */
	/*    o              */
	/* 	  4              */

	stree.insert_simplex_and_subfaces({0,1});
	stree.insert_simplex_and_subfaces({1,2,3});
	stree.insert_simplex_and_subfaces({2,4});
	stree.insert_simplex_and_subfaces({3,4});
	// Toy Simplicial Complex formed

	// Creating the MsMatrix named 'mat' from the Simplex_tree 'stree'
 	MsMatrix mat(stree);

 	// .collapsed_tree() does the strong collapse and returns the Simplex_tree of the collapsed Simplicial Complex
 	Simplex_tree stree_of_collapsed = mat.collapsed_tree();

 	// .reduction_map() returns the Reduction Map of Collapses
 	Map red_map = mat.reduction_map();
	
	std::cout << "Reduction Map is: " << std::endl;
	for (Map::iterator iter = red_map.begin(); iter != red_map.end(); iter++)
	{
		std::cout << "Vertex " << iter->first << " collapses to Vertex " << iter->second << std::endl;
	}
}
