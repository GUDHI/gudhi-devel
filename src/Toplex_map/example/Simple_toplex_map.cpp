/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2017
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

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Fake_simplex_tree.h>

#include <iostream>
#include <utility>  // for pair
#include <vector>

using Toplex_map = Gudhi::Fake_simplex_tree;
using typeVectorVertex = std::vector< Toplex_map::Vertex_handle >;
using typePairSimplexBool = std::pair< Toplex_map::Simplex_handle, bool >;

int main(int argc, char * const argv[]) {

    // TEST OF INSERTION
    std::cout << "********************************************************************" << std::endl;
    std::cout << "EXAMPLE OF SIMPLE INSERTION" << std::endl;
    // Construct the Toplex_map
    Toplex_map t_map;

    /* Simplex to be inserted:  */
    /*    1                     */
    /*    o                     */
    /*   /X\                    */
    /*  o---o---o               */
    /*  2   0   3               */

    // ++ FIRST
    std::cout << "   * INSERT 0" << std::endl;
    typeVectorVertex firstSimplexVector = { 0 };
    typePairSimplexBool returnValue = t_map.insert_simplex_and_subfaces(firstSimplexVector, 0.1);

    if (returnValue.second == true) {
        std::cout << "   + 0 INSERTED" << std::endl;
    } else {
        std::cout << "   - 0 NOT INSERTED" << std::endl;
    }

    // ++ SECOND
    std::cout << "   * INSERT 1" << std::endl;
    typeVectorVertex secondSimplexVector = { 1 };
    returnValue = t_map.insert_simplex_and_subfaces(secondSimplexVector, 0.1);

    if (returnValue.second == true) {
        std::cout << "   + 1 INSERTED" << std::endl;
    } else {
        std::cout << "   - 1 NOT INSERTED" << std::endl;
    }

    // ++ THIRD
    std::cout << "   * INSERT (0,1)" << std::endl;
    typeVectorVertex thirdSimplexVector = { 0, 1 };
    returnValue =
            t_map.insert_simplex_and_subfaces(thirdSimplexVector, 0.2);

    if (returnValue.second == true) {
        std::cout << "   + (0,1) INSERTED" << std::endl;
    } else {
        std::cout << "   - (0,1) NOT INSERTED" << std::endl;
    }

    // ++ FOURTH
    std::cout << "   * INSERT 2" << std::endl;
    typeVectorVertex fourthSimplexVector = { 2 };
    returnValue =
            t_map.insert_simplex_and_subfaces(fourthSimplexVector, 0.1);

    if (returnValue.second == true) {
        std::cout << "   + 2 INSERTED" << std::endl;
    } else {
        std::cout << "   - 2 NOT INSERTED" << std::endl;
    }

    // ++ FIFTH
    std::cout << "   * INSERT (2,0)" << std::endl;
    typeVectorVertex fifthSimplexVector = { 2, 0 };
    returnValue =
            t_map.insert_simplex_and_subfaces(fifthSimplexVector, 0.2);

    if (returnValue.second == true) {
        std::cout << "   + (2,0) INSERTED" << std::endl;
    } else {
        std::cout << "   - (2,0) NOT INSERTED" << std::endl;
    }

    // ++ SIXTH
    std::cout << "   * INSERT (2,1)" << std::endl;
    typeVectorVertex sixthSimplexVector = { 2, 1 };
    returnValue =
            t_map.insert_simplex_and_subfaces(sixthSimplexVector, 0.2);

    if (returnValue.second == true) {
        std::cout << "   + (2,1) INSERTED" << std::endl;
    } else {
        std::cout << "   - (2,1) NOT INSERTED" << std::endl;
    }

    // ++ SEVENTH
    std::cout << "   * INSERT (2,1,0)" << std::endl;
    typeVectorVertex seventhSimplexVector = { 2, 1, 0 };
    returnValue =
            t_map.insert_simplex_and_subfaces(seventhSimplexVector, 0.3);

    if (returnValue.second == true) {
        std::cout << "   + (2,1,0) INSERTED" << std::endl;
    } else {
        std::cout << "   - (2,1,0) NOT INSERTED" << std::endl;
    }

    // ++ EIGHTH
    std::cout << "   * INSERT 3" << std::endl;
    typeVectorVertex eighthSimplexVector = { 3 };
    returnValue =
            t_map.insert_simplex_and_subfaces(eighthSimplexVector, 0.1);

    if (returnValue.second == true) {
        std::cout << "   + 3 INSERTED" << std::endl;
    } else {
        std::cout << "   - 3 NOT INSERTED" << std::endl;
    }

    // ++ NINETH
    std::cout << "   * INSERT (3,0)" << std::endl;
    typeVectorVertex ninethSimplexVector = { 3, 0 };
    returnValue =
            t_map.insert_simplex_and_subfaces(ninethSimplexVector, 0.2);

    if (returnValue.second == true) {
        std::cout << "   + (3,0) INSERTED" << std::endl;
    } else {
        std::cout << "   - (3,0) NOT INSERTED" << std::endl;
    }

    // ++ TENTH
    std::cout << "   * INSERT 0 (already inserted)" << std::endl;
    typeVectorVertex tenthSimplexVector = { 0 };
    // With a different filtration value
    returnValue = t_map.insert_simplex_and_subfaces(tenthSimplexVector, 0.4);

    if (returnValue.second == true) {
        std::cout << "   + 0 INSERTED" << std::endl;
    } else {
        std::cout << "   - 0 NOT INSERTED" << std::endl;
    }

    // ++ ELEVENTH
    std::cout << "   * INSERT (2,1,0) (already inserted)" << std::endl;
    typeVectorVertex eleventhSimplexVector = { 2, 1, 0 };
    returnValue =
            t_map.insert_simplex_and_subfaces(eleventhSimplexVector, 0.4);

    if (returnValue.second == true) {
        std::cout << "   + (2,1,0) INSERTED" << std::endl;
    } else {
        std::cout << "   - (2,1,0) NOT INSERTED" << std::endl;
    }

    // ++ GENERAL VARIABLE SET

    std::cout << "********************************************************************\n";
    // Display the Simplex_tree - Can not be done in the middle of 2 inserts
    std::cout << "* The complex contains " << t_map.num_vertices() << " vertices and " << t_map.num_simplices()
              << " simplices - dimension is " << t_map.dimension() << "\n";
    std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
    for (auto f_simplex : t_map.filtration_simplex_range()) {
        std::cout << "   " << "[" << t_map.filtration(f_simplex) << "] ";
        for (auto vertex : t_map.simplex_vertex_range(f_simplex))
            std::cout << "(" << vertex << ")";
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

    std::cout << std::endl << std::endl;

    std::cout << "Iterator on skeleton:" << std::endl;
    for (auto f_simplex : t_map.skeleton_simplex_range()) {
        std::cout << "   " << "[" << t_map.filtration(f_simplex) << "] ";
        for (auto vertex : t_map.simplex_vertex_range(f_simplex)) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
