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
typedef Simplex_tree<> typeST;

int main(int argc, char * const argv[]) {
  // TEST OF INSERTION
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF INSERTION" << std::endl;
  typeST st;

  // ++ FIRST
  std::cout << "   - INSERT (2,1,0)" << std::endl;
  typeVectorVertex SimplexVector1;
  SimplexVector1.push_back(2);
  SimplexVector1.push_back(1);
  SimplexVector1.push_back(0);
  st.insert_simplex_and_subfaces(SimplexVector1, 0.3);

  // ++ SECOND
  std::cout << "   - INSERT 3" << std::endl;
  typeVectorVertex SimplexVector2;
  SimplexVector2.push_back(3);
  st.insert_simplex_and_subfaces(SimplexVector2, 0.1);

  // ++ THIRD
  std::cout << "   - INSERT (0,3)" << std::endl;
  typeVectorVertex SimplexVector3;
  SimplexVector3.push_back(3);
  SimplexVector3.push_back(0);
  st.insert_simplex_and_subfaces(SimplexVector3, 0.2);

  // ++ FOURTH
  std::cout << "   - INSERT (1,0) (already inserted)" << std::endl;
  typeVectorVertex SimplexVector4;
  SimplexVector4.push_back(1);
  SimplexVector4.push_back(0);
  st.insert_simplex_and_subfaces(SimplexVector4, 0.2);

  // ++ FIFTH
  std::cout << "   - INSERT (3,4,5)" << std::endl;
  typeVectorVertex SimplexVector5;
  SimplexVector5.push_back(3);
  SimplexVector5.push_back(4);
  SimplexVector5.push_back(5);
  st.insert_simplex_and_subfaces(SimplexVector5, 0.3);

  // ++ SIXTH
  std::cout << "   - INSERT (0,1,6,7)" << std::endl;
  typeVectorVertex SimplexVector6;
  SimplexVector6.push_back(0);
  SimplexVector6.push_back(1);
  SimplexVector6.push_back(6);
  SimplexVector6.push_back(7);
  st.insert_simplex_and_subfaces(SimplexVector6, 0.4);

  /* Inserted simplex:         */
  /*    1   6                  */
  /*    o---o                  */
  /*   /X\7/      4            */
  /*  o---o---o---o            */
  /*  2   0   3\X/             */
  /*            o              */
  /*            5              */

  /* In other words:          */
  /*   A facet [2,1,0]        */
  /*   An edge [0,3]          */
  /*   A facet [3,4,5]        */
  /*   A cell  [0,1,6,7]      */
  /*   A cell  [4,5,8,9]      */
  /*   A facet [9,10,11]      */
  /*   An edge [11,6]         */
  /*   An edge [10,12,2]      */

  // ++ GENERAL VARIABLE SET
  st.set_filtration(0.4); // Max filtration value
  st.set_dimension(3); // Max dimension = 3 -> (0,1,6,7)

  std::cout << "The complex contains " << st.num_simplices() << " simplices - " << st.num_vertices() << " vertices " << std::endl;
  std::cout << "   - dimension " << st.dimension() << "   - filtration " << st.filtration() << std::endl;
  std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  std::cout << "**************************************************************" << std::endl;

  for( auto f_simplex : st.filtration_simplex_range() )
  {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for( auto vertex : st.simplex_vertex_range(f_simplex) )
    {
      std::cout << (int)vertex;
    }
  }

  return 0;
}
