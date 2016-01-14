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
#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/Landmark_choice_by_random_point.h>
#include <gudhi/Landmark_choice_by_furthest_point.h>


using namespace Gudhi;

typedef std::vector< Vertex_handle > typeVectorVertex;
typedef Witness_complex<Simplex_tree<>> WitnessComplex;
typedef std::vector<double> Point;
//typedef std::pair<typeVectorVertex, Filtration_value> typeSimplex;
//typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;

int main (int argc, char * const argv[])
{
  std::vector< typeVectorVertex > knn;
  std::vector< Point > points;
  // Add grid points as witnesses
  for (double i = 0; i < 10; i += 1.0)
    for (double j = 0; j < 10; j += 1.0)
      for (double k = 0; k < 10; k += 1.0)
        points.push_back(Point({i,j,k}));

  bool b_print_output = false;
  // First test: random choice
  Simplex_tree<> complex1;
  Landmark_choice_by_random_point lcrp(points, 100, knn);
  assert(!knn.empty());
  WitnessComplex witnessComplex1(knn, complex1, 100, 3);
  assert(witnessComplex1.is_witness_complex(knn, b_print_output));

  // Second test: furthest choice
  knn.clear();
  Simplex_tree<> complex2;
  Landmark_choice_by_furthest_point lcfp(points, 100, knn);
  WitnessComplex witnessComplex2(knn, complex2, 100, 3);
  assert(witnessComplex2.is_witness_complex(knn, b_print_output));
}
