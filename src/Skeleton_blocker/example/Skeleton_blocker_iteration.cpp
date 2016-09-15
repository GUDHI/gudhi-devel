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

#include <boost/timer/timer.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>


using namespace std;
using namespace Gudhi;
using namespace skeleton_blocker;

typedef Skeleton_blocker_complex<Skeleton_blocker_simple_traits> Complex;
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
  boost::timer::auto_cpu_timer t;

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
  return EXIT_SUCCESS;
}
