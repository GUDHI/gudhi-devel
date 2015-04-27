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
typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;

int main(int argc, char * const argv[]) {
  const Filtration_value FIRST_FILTRATION_VALUE = 0.1;
  const Filtration_value SECOND_FILTRATION_VALUE = 0.2;
  const Filtration_value THIRD_FILTRATION_VALUE = 0.3;
  const Filtration_value FOURTH_FILTRATION_VALUE = 0.4;
  Vertex_handle FIRST_VERTEX_HANDLE = (Vertex_handle) 0;
  Vertex_handle SECOND_VERTEX_HANDLE = (Vertex_handle) 1;
  Vertex_handle THIRD_VERTEX_HANDLE = (Vertex_handle) 2;
  Vertex_handle FOURTH_VERTEX_HANDLE = (Vertex_handle) 3;

  // TEST OF INSERTION
  std::cout << "********************************************************************" << std::endl;
  std::cout << "EXAMPLE OF SIMPLE INSERTION" << std::endl;
  //Construct the Simplex Tree
  Simplex_tree<> simplexTree;

  /* Simplex to be inserted:  */
  /*    1                     */
  /*    o                     */
  /*   /X\                    */
  /*  o---o---o               */
  /*  2   0   3               */

  // ++ FIRST
  std::cout << "   * INSERT 0" << std::endl;
  typeVectorVertex firstSimplexVector;
  firstSimplexVector.push_back(FIRST_VERTEX_HANDLE);
  typePairSimplexBool returnValue =
      simplexTree.insert_simplex(firstSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + 0 INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - 0 NOT INSERTED" << std::endl;
  }

  // ++ SECOND
  std::cout << "   * INSERT 1" << std::endl;
  typeVectorVertex secondSimplexVector;
  secondSimplexVector.push_back(SECOND_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(secondSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + 1 INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - 1 NOT INSERTED" << std::endl;
  }

  // ++ THIRD
  std::cout << "   * INSERT (0,1)" << std::endl;
  typeVectorVertex thirdSimplexVector;
  thirdSimplexVector.push_back(FIRST_VERTEX_HANDLE);
  thirdSimplexVector.push_back(SECOND_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(thirdSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + (0,1) INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - (0,1) NOT INSERTED" << std::endl;
  }

  // ++ FOURTH
  std::cout << "   * INSERT 2" << std::endl;
  typeVectorVertex fourthSimplexVector;
  fourthSimplexVector.push_back(THIRD_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(fourthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + 2 INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - 2 NOT INSERTED" << std::endl;
  }

  // ++ FIFTH
  std::cout << "   * INSERT (2,0)" << std::endl;
  typeVectorVertex fifthSimplexVector;
  fifthSimplexVector.push_back(THIRD_VERTEX_HANDLE);
  fifthSimplexVector.push_back(FIRST_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(fifthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + (2,0) INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - (2,0) NOT INSERTED" << std::endl;
  }

  // ++ SIXTH
  std::cout << "   * INSERT (2,1)" << std::endl;
  typeVectorVertex sixthSimplexVector;
  sixthSimplexVector.push_back(THIRD_VERTEX_HANDLE);
  sixthSimplexVector.push_back(SECOND_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(sixthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + (2,1) INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - (2,1) NOT INSERTED" << std::endl;
  }

  // ++ SEVENTH
  std::cout << "   * INSERT (2,1,0)" << std::endl;
  typeVectorVertex seventhSimplexVector;
  seventhSimplexVector.push_back(THIRD_VERTEX_HANDLE);
  seventhSimplexVector.push_back(SECOND_VERTEX_HANDLE);
  seventhSimplexVector.push_back(FIRST_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(seventhSimplexVector, Filtration_value(THIRD_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + (2,1,0) INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - (2,1,0) NOT INSERTED" << std::endl;
  }

  // ++ EIGHTH
  std::cout << "   * INSERT 3" << std::endl;
  typeVectorVertex eighthSimplexVector;
  eighthSimplexVector.push_back(FOURTH_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(eighthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + 3 INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - 3 NOT INSERTED" << std::endl;
  }

  // ++ NINETH
  std::cout << "   * INSERT (3,0)" << std::endl;
  typeVectorVertex ninethSimplexVector;
  ninethSimplexVector.push_back(FOURTH_VERTEX_HANDLE);
  ninethSimplexVector.push_back(FIRST_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(ninethSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + (3,0) INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - (3,0) NOT INSERTED" << std::endl;
  }

  // ++ TENTH
  std::cout << "   * INSERT 0 (already inserted)" << std::endl;
  typeVectorVertex tenthSimplexVector;
  tenthSimplexVector.push_back(FIRST_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(tenthSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE)); // With a different filtration value

  if (returnValue.second == true) {
    std::cout << "   + 0 INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
    std::cout << "   - 0 NOT INSERTED" << std::endl;
  }

  // ++ ELEVENTH
  std::cout << "   * INSERT (2,1,0) (already inserted)" << std::endl;
  typeVectorVertex eleventhSimplexVector;
  eleventhSimplexVector.push_back(THIRD_VERTEX_HANDLE);
  eleventhSimplexVector.push_back(SECOND_VERTEX_HANDLE);
  eleventhSimplexVector.push_back(FIRST_VERTEX_HANDLE);
  returnValue =
      simplexTree.insert_simplex(eleventhSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::cout << "   + (2,1,0) INSERTED" << std::endl;
    int nb_simplices = simplexTree.num_simplices() + 1;
    simplexTree.set_num_simplices(nb_simplices);
  } else {
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
  for (auto f_simplex : simplexTree.filtration_simplex_range()) {
    std::cout << "   " << "[" << simplexTree.filtration(f_simplex) << "] ";
    for (auto vertex : simplexTree.simplex_vertex_range(f_simplex)) {
      std::cout << (int) vertex << " ";
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

  // ------------------------------------------------------------------------------------------------------------------
  // Find in the simplex_tree
  // ------------------------------------------------------------------------------------------------------------------
  Simplex_tree<>::Simplex_handle simplexFound = simplexTree.find(secondSimplexVector);
  std::cout << "**************IS THE SIMPLEX {1} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";

  Vertex_handle UNKNOWN_VERTEX_HANDLE = (Vertex_handle) 15;
  typeVectorVertex unknownSimplexVector;
  unknownSimplexVector.push_back(UNKNOWN_VERTEX_HANDLE);
  simplexFound = simplexTree.find(unknownSimplexVector);
  std::cout << "**************IS THE SIMPLEX {15} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";

  simplexFound = simplexTree.find(fifthSimplexVector);
  std::cout << "**************IS THE SIMPLEX {2,0} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";

  typeVectorVertex otherSimplexVector;
  otherSimplexVector.push_back(UNKNOWN_VERTEX_HANDLE);
  otherSimplexVector.push_back(SECOND_VERTEX_HANDLE);
  simplexFound = simplexTree.find(otherSimplexVector);
  std::cout << "**************IS THE SIMPLEX {15,1} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";

  typeVectorVertex invSimplexVector;
  invSimplexVector.push_back(SECOND_VERTEX_HANDLE);
  invSimplexVector.push_back(THIRD_VERTEX_HANDLE);
  invSimplexVector.push_back(FIRST_VERTEX_HANDLE);
  simplexFound = simplexTree.find(invSimplexVector);
  std::cout << "**************IS THE SIMPLEX {1,2,0} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";
  return 0;
}
