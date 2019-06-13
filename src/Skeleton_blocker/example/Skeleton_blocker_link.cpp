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

typedef Gudhi::skeleton_blocker::Skeleton_blocker_simple_traits Traits;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_complex<Traits> Complex;
typedef Complex::Vertex_handle Vertex_handle;
typedef Complex::Root_vertex_handle Root_vertex_handle;
typedef Complex::Simplex Simplex;

int main(int argc, char *argv[]) {
  // build a full complex with 4 vertices and 2^4-1 simplices

  // Create a complex with four vertices (0,1,2,3)
  Complex complex;

  // Add a tetrahedron to this complex
  Simplex tetrahedron(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2), Vertex_handle(3));
  complex.add_simplex(tetrahedron);

  std::cout << "complex:" << complex.to_string() << std::endl;

  // build the link of vertex 1, eg a triangle {0,2,3}
  auto link = complex.link(Vertex_handle(1));
  std::cout << "link:" << link.to_string() << std::endl;

  // Internally link is a subcomplex of 'complex' and its vertices are stored in a vector.
  // They can be accessed via Vertex_handle(x) where x is an index of the vector.
  // In that example, link has three vertices and thus it contains only
  // Vertex_handle(0),Vertex_handle(1) and Vertex_handle(2) are).
  for (int i = 0; i < 5; ++i)
    std::cout << "link.contains_vertex(Vertex_handle(" << i << ")):" << link.contains_vertex(Vertex_handle(i)) <<
        std::endl;
  std::cout << std::endl;

  // To access to the initial vertices eg (0,1,2,3,4),  Root_vertex_handle must be used.
  // For instance, to test if the link contains the vertex that was labeled i:
  for (int i = 0; i < 5; ++i)
    std::cout << "link.contains_vertex(Root_vertex_handle(" << i << ")):" <<
      link.contains_vertex(Root_vertex_handle(i)) << std::endl;

  return EXIT_SUCCESS;
}
