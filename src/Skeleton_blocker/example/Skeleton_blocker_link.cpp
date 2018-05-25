/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
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
