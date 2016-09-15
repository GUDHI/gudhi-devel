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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include "gudhi/Test.h"
#include "gudhi/Skeleton_blocker.h"


using namespace std;
using namespace Gudhi;
using namespace skeleton_blocker;

struct Geometry_trait {
  typedef std::vector<double> Point;
};

typedef Geometry_trait::Point Point;
typedef Skeleton_blocker_simple_geometric_traits<Geometry_trait> Complex_geometric_traits;
typedef Skeleton_blocker_geometric_complex< Complex_geometric_traits > Complex;
typedef Complex::Vertex_handle Vertex_handle;

bool test_constructor1() {
  Complex complex;
  Skeleton_blocker_off_reader<Complex> off_reader("test2.off", complex);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file" << std::endl;
    return false;
  }


  std::cout << "complex has " <<
      complex.num_vertices() << " vertices, " <<
      complex.num_blockers() << " blockers, " <<
      complex.num_edges() << " edges and" <<
      complex.num_triangles() << " triangles.";

  if (complex.num_vertices() != 7 || complex.num_edges() != 12 || complex.num_triangles() != 6)
    return false;

  Skeleton_blocker_off_writer<Complex> off_writer("tmp.off", complex);
  Complex same;
  Skeleton_blocker_off_reader<Complex> off_reader2("tmp.off", same);

  std::cout << "\ncomplex:" << complex.to_string() << endl;
  std::cout << "\nsame:" << same.to_string() << endl;

  return (complex == same);
}

bool test_constructor2() {
  Complex complex;
  Skeleton_blocker_off_reader<Complex> off_reader("test2.off", complex);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file" << std::endl;
    return false;
  }
  std::cout << "complex has " <<
      complex.num_vertices() << " vertices, " <<
      complex.num_blockers() << " blockers, " <<
      complex.num_edges() << " edges and" <<
      complex.num_triangles() << " triangles.";

  if (complex.num_vertices() != 7 || complex.num_edges() != 12 || complex.num_triangles() != 6)
    return false;

  auto link_0 = complex.abstract_link(Vertex_handle(0));


  std::cout << "\n link(0):" << link_0.to_string() << endl;

  auto link_geometric_0 = complex.link(Vertex_handle(0));

  auto print_point = [&](Vertex_handle v) {
    for (auto x : link_geometric_0.point(v)) std::cout << x << " ";
    std::cout << std::endl;
  };

  std::for_each(link_geometric_0.vertex_range().begin(), link_geometric_0.vertex_range().end(), print_point);

  //	for(auto v : link_geometric_0.vertex_range())
  //		std::cout<<"point("<<v<<"):"<<link_geometric_0.point(v)<<std::endl;

  return link_0.num_vertices() == 2;
}

int main(int argc, char *argv[]) {
  Tests tests_geometric_complex;
  tests_geometric_complex.add("Test constructor 1", test_constructor1);
  tests_geometric_complex.add("Test constructor 2", test_constructor2);

  if (tests_geometric_complex.run())
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
