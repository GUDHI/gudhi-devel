/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA
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

#include <gudhi/Skeleton_blocker.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace Gudhi;
using namespace skeleton_blocker;

typedef Skeleton_blocker_complex<Skeleton_blocker_simple_traits> Complex;
typedef Complex::Vertex_handle Vertex_handle;
typedef Complex::Simplex Simplex;

int main(int argc, char *argv[]) {
  std::vector<Simplex> simplices;

  // add 4 triangles of a tetrahedron 0123
  simplices.push_back(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2)));
  simplices.push_back(Simplex(Vertex_handle(1), Vertex_handle(2), Vertex_handle(3)));
  simplices.push_back(Simplex(Vertex_handle(3), Vertex_handle(0), Vertex_handle(2)));
  simplices.push_back(Simplex(Vertex_handle(3), Vertex_handle(0), Vertex_handle(1)));

  // get complex from top faces
  Complex complex(make_complex_from_top_faces<Complex>(simplices.begin(), simplices.end()));


  std::cout << "Simplices:" << std::endl;
  for (const Simplex & s : complex.complex_simplex_range())
    std::cout << s << " ";
  std::cout << std::endl;

  // One blocker as simplex 0123 is not in the complex but all its proper faces are.
  std::cout << "Blockers: " << complex.blockers_to_string() << std::endl;

  // now build a complex from its full list of simplices
  simplices.clear();
  simplices.push_back(Simplex(Vertex_handle(0)));
  simplices.push_back(Simplex(Vertex_handle(1)));
  simplices.push_back(Simplex(Vertex_handle(2)));
  simplices.push_back(Simplex(Vertex_handle(0), Vertex_handle(1)));
  simplices.push_back(Simplex(Vertex_handle(1), Vertex_handle(2)));
  simplices.push_back(Simplex(Vertex_handle(2), Vertex_handle(0)));
  complex = Complex(simplices.begin(), simplices.end());

  std::cout << "Simplices:" << std::endl;
  for (const Simplex & s : complex.complex_simplex_range())
    std::cout << s << " ";
  std::cout << std::endl;

  // One blocker as simplex 012 is not in the complex but all its proper faces are.
  std::cout << "Blockers: " << complex.blockers_to_string() << std::endl;

  return EXIT_SUCCESS;
}
