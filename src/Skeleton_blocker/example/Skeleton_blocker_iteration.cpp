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
#include <gudhi/Clock.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>

typedef Gudhi::skeleton_blocker::Skeleton_blocker_simple_traits Traits;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_complex<Traits> Complex;
typedef Complex::Vertex_handle Vertex_handle;
typedef Complex::Simplex Simplex;

Complex build_complete_complex(int n) {
  // build a full complex with n vertices and 2^n-1 simplices
  Complex complex;
  for (int i = 0; i < n; i++)
    complex.add_vertex();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
  return complex;
}

int main(int argc, char *argv[]) {
  Gudhi::Clock skbl_chrono("Time to build the complete complex, enumerate simplices and Euler Characteristic");
  const int n = 15;

  // build a full complex with n vertices and 2^n-1 simplices
  Complex complex(build_complete_complex(n));

  // this is just to illustrate iterators, to count number of vertices
  // or edges, complex.num_vertices() and complex.num_edges() are
  // more appropriated!
  unsigned num_vertices = 0;
  for (auto v : complex.vertex_range()) {
    std::cout << "Vertex " << v << std::endl;
    ++num_vertices;
  }

  // such loop can also be done directly with distance as iterators are STL compliant
  auto edges = complex.edge_range();
  unsigned num_edges = std::distance(edges.begin(), edges.end());

  unsigned euler = 0;
  unsigned num_simplices = 0;
  // we use a reference to a simplex instead of a copy
  // value here because a simplex is a set of integers
  // and copying it cost time
  for (const Simplex & s : complex.complex_simplex_range()) {
    ++num_simplices;
    if (s.dimension() % 2 == 0)
      euler += 1;
    else
      euler -= 1;
  }
  std::cout << "Saw " << num_vertices << " vertices, " << num_edges << " edges and " << num_simplices << " simplices"
      << std::endl;
  std::cout << "The Euler Characteristic is " << euler << std::endl;
  std::cout << skbl_chrono;
  return EXIT_SUCCESS;
}
