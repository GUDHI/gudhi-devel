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

	int tot;
	std::cout << "Enter the number of MaxSimplices in the Simplicial Complex you want to collapse :" << std::endl;
	std::cin >> tot;

	std::cout << "Now, enter those MaxSimplices line by line. The format of one MaxSimplex is :" << std::endl;
	std::cout << "k v_1 .... v_k" << std::endl;
	std::cout << "k represents number of vertices in this Simplex. Follow it with the k vertices of this MaxSimplex." << std::endl;

	int num_vert_in_this_mx_simplex;
	int vertex;
	typeVectorVertex to_ins;

	for(int run = 0 ; run < tot ; ++run)
	{
		std::cin >> num_vert_in_this_mx_simplex;
		to_ins.clear();
		for(int in = 0 ; in < num_vert_in_this_mx_simplex ; ++in)
		{
			std::cin >> vertex;
			to_ins.push_back(vertex);
		}
		stree.insert_simplex_and_subfaces(to_ins);
 	}

 	clock_t stree_formed = clock();
 	std::cout << "Simplex tree formed ... Now going for collapse" << std::endl;

 	MsMatrix mat(stree);
 	clock_t matrix_formed = clock();

 	// .strong_collapse() does strong collapse without computing the Simplex Tree corresponding to the Resultant Simplicial Complex
 	mat.strong_collapse();
 	clock_t collapse_done = clock();

 	std::cout << "Collapse done !" << std::endl;
 	
 	Map red_map = mat.reduction_map();
	
	std::cout << "Reduction Map is: " << std::endl;
	for (Map::iterator iter = red_map.begin(); iter != red_map.end(); iter++)
	{
		std::cout << "Key: " << iter->first << " and " << "Value: " << iter->second << std::endl;
	}

	std::cout << "Time for formation of Matrix : " << (matrix_formed - stree_formed)/CLOCKS_PER_SEC  << " seconds." << std::endl;
	std::cout << "Time for Collapse : " << (collapse_done - matrix_formed)/CLOCKS_PER_SEC  << " seconds." << std::endl;
}
