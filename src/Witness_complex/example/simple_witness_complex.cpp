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
//#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Simplex_tree.h"
#include "gudhi/Witness_complex.h"

using namespace Gudhi;

typedef std::vector< Vertex_handle > typeVectorVertex;
typedef Witness_complex<Simplex_tree<>> WitnessComplex;
//typedef std::pair<typeVectorVertex, Filtration_value> typeSimplex;
//typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;

int main (int argc, char * const argv[])
{
  Simplex_tree<> complex;
  std::vector< typeVectorVertex > knn;
  typeVectorVertex witness0  = {1,0,5,2,6,3,4}; knn.push_back(witness0 );
  typeVectorVertex witness1  = {2,6,4,5,0,1,3}; knn.push_back(witness1 );
  typeVectorVertex witness2  = {3,4,2,1,5,6,0}; knn.push_back(witness2 );
  typeVectorVertex witness3  = {4,2,1,3,5,6,0}; knn.push_back(witness3 );
  typeVectorVertex witness4  = {5,1,6,0,2,3,4}; knn.push_back(witness4 );
  typeVectorVertex witness5  = {6,0,5,2,1,3,4}; knn.push_back(witness5 );
  typeVectorVertex witness6  = {0,5,6,1,2,3,4}; knn.push_back(witness6 );
  typeVectorVertex witness7  = {2,6,4,5,3,1,0}; knn.push_back(witness7 );
  typeVectorVertex witness8  = {1,2,5,4,3,6,0}; knn.push_back(witness8 );
  typeVectorVertex witness9  = {3,4,0,6,5,1,2}; knn.push_back(witness9 );
  typeVectorVertex witness10 = {5,0,1,3,6,2,4}; knn.push_back(witness10);
  typeVectorVertex witness11 = {5,6,1,0,2,3,4}; knn.push_back(witness11);
  typeVectorVertex witness12 = {1,6,0,5,2,3,4}; knn.push_back(witness12);
  WitnessComplex witnessComplex(knn, complex, 7, 7);
  if (witnessComplex.is_witness_complex(knn, true))
    std::cout << "Witness complex is good\n";
  else
    std::cout << "Witness complex is bad\n";
}
