/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <utility>  // for pair
#include <vector>

using Simplex_tree = Gudhi::Simplex_tree<>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using typeVectorVertex = std::vector<Vertex_handle>;
using typePairSimplexBool = std::pair<Simplex_tree::Simplex_handle, bool>;

int main(int argc, char* const argv[]) {
  const Filtration_value FIRST_FILTRATION_VALUE = 0.1;
  const Filtration_value SECOND_FILTRATION_VALUE = 0.2;
  const Filtration_value THIRD_FILTRATION_VALUE = 0.3;
  const Filtration_value FOURTH_FILTRATION_VALUE = 0.4;

  // TEST OF INSERTION
  std::clog << "********************************************************************" << std::endl;
  std::clog << "EXAMPLE OF SIMPLE INSERTION" << std::endl;
  // Construct the Simplex Tree
  Simplex_tree simplexTree;

  /* Simplex to be inserted:  */
  /*    1                     */
  /*    o                     */
  /*   /X\                    */
  /*  o---o---o               */
  /*  2   0   3               */

  // ++ FIRST
  std::clog << "   * INSERT 0" << std::endl;
  typeVectorVertex firstSimplexVector = {0};
  typePairSimplexBool returnValue =
      simplexTree.insert_simplex(firstSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + 0 INSERTED" << std::endl;
  } else {
    std::clog << "   - 0 NOT INSERTED" << std::endl;
  }

  // ++ SECOND
  std::clog << "   * INSERT 1" << std::endl;
  typeVectorVertex secondSimplexVector = {1};
  returnValue = simplexTree.insert_simplex(secondSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + 1 INSERTED" << std::endl;
  } else {
    std::clog << "   - 1 NOT INSERTED" << std::endl;
  }

  // ++ THIRD
  std::clog << "   * INSERT (0,1)" << std::endl;
  typeVectorVertex thirdSimplexVector = {0, 1};
  returnValue = simplexTree.insert_simplex(thirdSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + (0,1) INSERTED" << std::endl;
  } else {
    std::clog << "   - (0,1) NOT INSERTED" << std::endl;
  }

  // ++ FOURTH
  std::clog << "   * INSERT 2" << std::endl;
  typeVectorVertex fourthSimplexVector = {2};
  returnValue = simplexTree.insert_simplex(fourthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + 2 INSERTED" << std::endl;
  } else {
    std::clog << "   - 2 NOT INSERTED" << std::endl;
  }

  // ++ FIFTH
  std::clog << "   * INSERT (2,0)" << std::endl;
  typeVectorVertex fifthSimplexVector = {2, 0};
  returnValue = simplexTree.insert_simplex(fifthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + (2,0) INSERTED" << std::endl;
  } else {
    std::clog << "   - (2,0) NOT INSERTED" << std::endl;
  }

  // ++ SIXTH
  std::clog << "   * INSERT (2,1)" << std::endl;
  typeVectorVertex sixthSimplexVector = {2, 1};
  returnValue = simplexTree.insert_simplex(sixthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + (2,1) INSERTED" << std::endl;
  } else {
    std::clog << "   - (2,1) NOT INSERTED" << std::endl;
  }

  // ++ SEVENTH
  std::clog << "   * INSERT (2,1,0)" << std::endl;
  typeVectorVertex seventhSimplexVector = {2, 1, 0};
  returnValue = simplexTree.insert_simplex(seventhSimplexVector, Filtration_value(THIRD_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + (2,1,0) INSERTED" << std::endl;
  } else {
    std::clog << "   - (2,1,0) NOT INSERTED" << std::endl;
  }

  // ++ EIGHTH
  std::clog << "   * INSERT 3" << std::endl;
  typeVectorVertex eighthSimplexVector = {3};
  returnValue = simplexTree.insert_simplex(eighthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + 3 INSERTED" << std::endl;
  } else {
    std::clog << "   - 3 NOT INSERTED" << std::endl;
  }

  // ++ NINETH
  std::clog << "   * INSERT (3,0)" << std::endl;
  typeVectorVertex ninethSimplexVector = {3, 0};
  returnValue = simplexTree.insert_simplex(ninethSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + (3,0) INSERTED" << std::endl;
  } else {
    std::clog << "   - (3,0) NOT INSERTED" << std::endl;
  }

  // ++ TENTH
  std::clog << "   * INSERT 0 (already inserted)" << std::endl;
  typeVectorVertex tenthSimplexVector = {0};
  // With a different filtration value
  returnValue = simplexTree.insert_simplex(tenthSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + 0 INSERTED" << std::endl;
  } else {
    std::clog << "   - 0 NOT INSERTED" << std::endl;
  }

  // ++ ELEVENTH
  std::clog << "   * INSERT (2,1,0) (already inserted)" << std::endl;
  typeVectorVertex eleventhSimplexVector = {2, 1, 0};
  returnValue = simplexTree.insert_simplex(eleventhSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE));

  if (returnValue.second == true) {
    std::clog << "   + (2,1,0) INSERTED" << std::endl;
  } else {
    std::clog << "   - (2,1,0) NOT INSERTED" << std::endl;
  }

  // ++ GENERAL VARIABLE SET

  std::clog << "********************************************************************\n";
  std::clog << "* The complex contains " << simplexTree.num_simplices() << " simplices\n";
  std::clog << "   - dimension " << simplexTree.dimension() << "\n";
  std::clog << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplexTree.filtration_simplex_range()) {
    std::clog << "   "
              << "[" << simplexTree.filtration(f_simplex) << "] ";
    for (auto vertex : simplexTree.simplex_vertex_range(f_simplex)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
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
  Simplex_tree::Simplex_handle simplexFound = simplexTree.find(secondSimplexVector);
  std::clog << "**************IS THE SIMPLEX {1} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";

  typeVectorVertex unknownSimplexVector = {15};
  simplexFound = simplexTree.find(unknownSimplexVector);
  std::clog << "**************IS THE SIMPLEX {15} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";

  simplexFound = simplexTree.find(fifthSimplexVector);
  std::clog << "**************IS THE SIMPLEX {2,0} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";

  typeVectorVertex otherSimplexVector = {1, 15};
  simplexFound = simplexTree.find(otherSimplexVector);
  std::clog << "**************IS THE SIMPLEX {15,1} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";

  typeVectorVertex invSimplexVector = {1, 2, 0};
  simplexFound = simplexTree.find(invSimplexVector);
  std::clog << "**************IS THE SIMPLEX {1,2,0} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";

  simplexFound = simplexTree.find({0, 1});
  std::clog << "**************IS THE SIMPLEX {0,1} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != simplexTree.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";

  std::clog << "**************COFACES OF {0,1} IN CODIMENSION 1 ARE\n";
  for (auto& simplex : simplexTree.cofaces_simplex_range(simplexTree.find({0, 1}), 1)) {
    for (auto vertex : simplexTree.simplex_vertex_range(simplex)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

  std::clog << "**************STARS OF {0,1} ARE\n";
  for (auto& simplex : simplexTree.star_simplex_range(simplexTree.find({0, 1}))) {
    for (auto vertex : simplexTree.simplex_vertex_range(simplex)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

  std::clog << "**************BOUNDARIES OF {0,1,2} ARE\n";
  for (auto& simplex : simplexTree.boundary_simplex_range(simplexTree.find({0, 1, 2}))) {
    for (auto vertex : simplexTree.simplex_vertex_range(simplex)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }

  return 0;
}
