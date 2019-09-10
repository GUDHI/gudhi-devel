/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Skeleton_blocker.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

typedef Gudhi::skeleton_blocker::Skeleton_blocker_simple_traits Traits;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_complex<Traits> Complex;
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
  Complex complex(Gudhi::skeleton_blocker::make_complex_from_top_faces<Complex>(simplices.begin(), simplices.end()));


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
