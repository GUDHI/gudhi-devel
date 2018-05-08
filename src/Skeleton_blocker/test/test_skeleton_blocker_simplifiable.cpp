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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "skeleton_blocker_simplifiable"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Skeleton_blocker.h>

template<typename ComplexType> class Skeleton_blocker_sub_complex;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_simple_traits Traits;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_complex<Traits> Complex;
typedef Complex::Vertex_handle Vertex_handle;
typedef Complex::Root_vertex_handle Root_vertex_handle;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_simplex<Vertex_handle> Simplex;


void build_complete(int n, Complex& complex) {
  complex.clear();
  for (int i = 0; i < n; i++)
    complex.add_vertex();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_contraction1) {
  enum { a, b, x, y, z, n };
  Complex complex(n);
  build_complete(n, complex);
  complex.remove_edge(static_cast<Vertex_handle> (b), static_cast<Vertex_handle> (z));
  complex.add_blocker(Simplex(static_cast<Vertex_handle> (a), static_cast<Vertex_handle> (x),
                              static_cast<Vertex_handle> (y)));
  complex.add_blocker(Simplex(static_cast<Vertex_handle> (b), static_cast<Vertex_handle> (x),
                              static_cast<Vertex_handle> (y)));

  // Print result
  std::cout << "complex before complex" << complex.to_string() << std::endl;

  std::cout << std::endl << std::endl;
  complex.contract_edge(static_cast<Vertex_handle> (a), static_cast<Vertex_handle> (b));
  // Print result
  std::cout << "ContractEdge(0,1)\n";
  PRINT(complex.to_string());

  // verification
  for (int i = 0; i < 5; i++)
    if (i != 1) BOOST_CHECK(complex.contains(Simplex(static_cast<Vertex_handle> (i))));
  BOOST_CHECK(!complex.contains_edge(static_cast<Vertex_handle> (a), static_cast<Vertex_handle> (b)));
  
  BOOST_CHECK(complex.contains_blocker(Simplex(*complex.get_address(Root_vertex_handle(a)),
                                               *complex.get_address(Root_vertex_handle(x)),
                                               *complex.get_address(Root_vertex_handle(y)))));
  
  BOOST_CHECK(complex.num_edges() == 6);
  BOOST_CHECK(complex.num_blockers() == 1);
  Simplex sigma;
  sigma.add_vertex(static_cast<Vertex_handle> (a));
  sigma.add_vertex(static_cast<Vertex_handle> (x));
  sigma.add_vertex(static_cast<Vertex_handle> (y));
  sigma.add_vertex(static_cast<Vertex_handle> (z));
  BOOST_CHECK(!(complex.contains(sigma)));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_contraction2) {
  enum { a, b, x, y, z, n };
  Complex complex(n);
  build_complete(n, complex);
  complex.remove_edge(static_cast<Vertex_handle> (b), static_cast<Vertex_handle> (x));
  Simplex blocker;
  blocker.add_vertex(static_cast<Vertex_handle> (a));
  blocker.add_vertex(static_cast<Vertex_handle> (y));
  blocker.add_vertex(static_cast<Vertex_handle> (z));

  complex.add_blocker(blocker);

  // Print result
  std::cout << "complex complex" << complex.to_string();
  std::cout << std::endl << std::endl;
  complex.contract_edge(static_cast<Vertex_handle> (a), static_cast<Vertex_handle> (b));

  std::cout << "complex.ContractEdge(a,b)" << complex.to_string();

  std::cout << std::endl << std::endl;

  // there should be one blocker (a,c,d,e) in the complex
  BOOST_CHECK(complex.contains_blocker(Simplex(static_cast<Vertex_handle> (a), static_cast<Vertex_handle> (x),
                                               static_cast<Vertex_handle> (y), static_cast<Vertex_handle> (z))));
  BOOST_CHECK(complex.num_blockers() == 1);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_link_condition1) {
  Complex complex(0);
  // Build the complexes
  build_complete(4, complex);
  complex.add_blocker(Simplex(static_cast<Vertex_handle> (0), static_cast<Vertex_handle> (1), static_cast<Vertex_handle> (2)));

  // Print result
  std::cout << "complex complex" << complex.to_string();
  std::cout << std::endl << std::endl;

  BOOST_CHECK(complex.link_condition(Vertex_handle(1), Vertex_handle(2), true));

  BOOST_CHECK(!complex.link_condition(Vertex_handle(1), Vertex_handle(2), false));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_collapse0) {
  Complex complex(5);
  build_complete(4, complex);
  complex.add_vertex();
  complex.add_edge_without_blockers(static_cast<Vertex_handle> (2), static_cast<Vertex_handle> (4));
  complex.add_edge_without_blockers(static_cast<Vertex_handle> (3), static_cast<Vertex_handle> (4));
  // Print result
  std::cout << "initial complex :\n" << complex.to_string();
  std::cout << std::endl << std::endl;

  Simplex simplex_123(static_cast<Vertex_handle> (1), static_cast<Vertex_handle> (2), static_cast<Vertex_handle> (3));
  complex.remove_star(simplex_123);
  std::cout << "complex.remove_star(1,2,3):\n" << complex.to_string();
  std::cout << std::endl << std::endl;

  // verification
  BOOST_CHECK(complex.contains_blocker(simplex_123));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_collapse1) {
  Complex complex(5);
  build_complete(4, complex);
  complex.add_blocker(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2), Vertex_handle(3)));
  // Print result
  std::cout << "initial complex :\n" << complex.to_string();
  std::cout << std::endl << std::endl;

  Simplex simplex_123(Vertex_handle(1), Vertex_handle(2), Vertex_handle(3));
  complex.remove_star(simplex_123);
  std::cout << "complex.remove_star(1,2,3):\n" << complex.to_string();
  std::cout << std::endl << std::endl;

  // verification
  BOOST_CHECK(complex.contains_blocker(simplex_123));
  BOOST_CHECK(complex.num_blockers() == 1);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_collapse2) {
  Complex complex(5);
  build_complete(4, complex);
  complex.add_vertex();
  complex.add_edge_without_blockers(Vertex_handle(1), Vertex_handle(4));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(4));
  complex.add_edge_without_blockers(Vertex_handle(3), Vertex_handle(4));
  complex.add_blocker(Simplex(Vertex_handle(1), Vertex_handle(2), Vertex_handle(3), Vertex_handle(4)));
  // Print result
  std::cout << "initial complex :\n" << complex.to_string();
  std::cout << std::endl << std::endl;

  Simplex sigma(Vertex_handle(1), Vertex_handle(2), Vertex_handle(3));
  complex.remove_star(sigma);
  std::cout << "complex.remove_star(1,2,3):\n" << complex.to_string();
  std::cout << std::endl << std::endl;

  // verification
  BOOST_CHECK(!complex.contains_blocker(Simplex(Vertex_handle(1), Vertex_handle(2),
                                                Vertex_handle(3), Vertex_handle(4))));
  BOOST_CHECK(complex.contains_blocker(sigma));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_collapse3) {
  Complex complex(5);
  build_complete(4, complex);
  complex.add_vertex();
  complex.add_edge_without_blockers(Vertex_handle(1), Vertex_handle(4));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(4));
  complex.add_edge_without_blockers(Vertex_handle(3), Vertex_handle(4));
  complex.add_blocker(Simplex(Vertex_handle(1), Vertex_handle(2), Vertex_handle(3), Vertex_handle(4)));
  // Print result
  std::cout << "initial complex:\n" << complex.to_string();
  std::cout << std::endl << std::endl;

  complex.remove_star(static_cast<Vertex_handle> (2));
  std::cout << "complex after remove star of 2:\n" << complex.to_string();

  BOOST_CHECK(complex.contains_blocker(Simplex(Vertex_handle(1), Vertex_handle(3), Vertex_handle(4))));
  BOOST_CHECK(!complex.contains_blocker(Simplex(Vertex_handle(1), Vertex_handle(2),
                                                Vertex_handle(3), Vertex_handle(4))));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_add_simplex) {
  Complex complex(4);
  build_complete(4, complex);
  complex.add_blocker(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(3)));
  std::cout << "initial complex:\n" << complex.to_string();
  std::cout << std::endl << std::endl;

  complex.add_simplex(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(3)));
  std::cout << "complex after add_simplex:\n" << complex.to_string();
  BOOST_CHECK(complex.num_blockers() == 1);
  BOOST_CHECK(complex.contains_blocker(Simplex(Vertex_handle(0), Vertex_handle(1),
                                               Vertex_handle(2), Vertex_handle(3))));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_add_simplex2) {
  Complex complex;
  build_complete(4, complex);
  // Print result
  std::cout << "initial complex:\n" << complex.to_string();
  std::cout << std::endl << std::endl;

  Complex copy(complex.num_vertices());

  std::vector<Simplex> simplices(complex.complex_simplex_range().begin(), complex.complex_simplex_range().end());
  sort(simplices.begin(), simplices.end(), [&](const Simplex& s1, const Simplex & s2) {
    return s1.dimension() < s2.dimension();
  });
  for (const auto & simplex : simplices) {
    if (!copy.contains(simplex) && simplex.dimension() == 1)
      copy.add_edge_without_blockers(simplex.first_vertex(), simplex.last_vertex());
    if (!copy.contains(simplex) && simplex.dimension() > 1)
      copy.add_simplex(simplex);
  }

  std::cout << "complex after add_simplex:\n" << copy.to_string();

  BOOST_CHECK(complex.num_blockers() == copy.num_blockers());
  BOOST_CHECK(complex.num_edges() == copy.num_edges());
  BOOST_CHECK(complex.num_vertices() == copy.num_vertices());
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_add_simplex3) {
  Complex complex(5);
  build_complete(5, complex);
  complex.remove_edge(Vertex_handle(3), Vertex_handle(4));
  Simplex sigma(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2));
  complex.add_blocker(sigma);
  // Print result
  std::cout << "initial complex:\n" << complex.to_string();
  std::cout << std::endl << std::endl;
  complex.add_simplex(sigma);
  //should create two blockers 0123 and 0124
  std::cout << "complex after adding simplex 012:\n" << complex.to_string();
  BOOST_CHECK(complex.num_blockers() == 2);
  BOOST_CHECK(complex.contains_blocker(Simplex(Vertex_handle(0), Vertex_handle(1),
                                               Vertex_handle(2), Vertex_handle(3))));
  BOOST_CHECK(complex.contains_blocker(Simplex(Vertex_handle(0), Vertex_handle(1),
                                               Vertex_handle(2), Vertex_handle(4))));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_add_simplex4) {
  int n = 6;
  Complex complex(n);

  // add all simplex 0..n without i
  for (int i = 0; i < n; i++) {
    Simplex s;
    for (int k = 0; k < n; k++)
      s.add_vertex(Vertex_handle(k));
    s.remove_vertex(Vertex_handle(i));
    complex.add_simplex(s);

    //at step i there is only blocker 0..i
    BOOST_CHECK(!(i < 2 && complex.num_blockers() > 0));
    if (i >= 2 && complex.num_blockers() != 1) {
      Simplex b;
      for (int k = 0; k < i; k++)
        b.add_vertex(Vertex_handle(i));
      BOOST_CHECK(complex.contains_blocker(b));
    }
  }
  Simplex s;
  for (int k = 0; k < n; k++)
    s.add_vertex(Vertex_handle(k));
  BOOST_CHECK(complex.num_blockers() == 1);
  BOOST_CHECK(complex.contains_blocker(s));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_add_edge) {
  Complex complex(4);
  for (unsigned i = 0u; i < 4; i++)
    complex.add_edge(Vertex_handle(i), Vertex_handle((i + 1) % 4));

  // Print result
  std::cout << "initial complex:\n" << complex.to_string();
  std::cout << std::endl << std::endl;
  complex.add_edge(Vertex_handle(1), Vertex_handle(3));
  //should create two blockers 013 and 012
  std::cout << "complex after adding edge 13:\n" << complex.to_string();
  BOOST_CHECK(complex.num_blockers() == 2);
  BOOST_CHECK(complex.contains_blocker(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(3))));
  BOOST_CHECK(complex.contains_blocker(Simplex(Vertex_handle(1), Vertex_handle(2), Vertex_handle(3))));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplifiable_remove_popable_blockers) {
  Complex complex;
  build_complete(4, complex);
  complex.add_vertex();
  complex.add_edge_without_blockers(Vertex_handle(3), Vertex_handle(4));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(4));
  Simplex sigma1 = Simplex(Vertex_handle(1), Vertex_handle(2), Vertex_handle(3));
  Simplex sigma2 = Simplex(Vertex_handle(2), Vertex_handle(3), Vertex_handle(4));

  complex.add_blocker(sigma1);
  complex.add_blocker(sigma2);
  std::cout << "complex complex" << complex.to_string();
  std::cout << std::endl << std::endl;
  std::cout << "complex.RemovePopableBlockers();" << std::endl;
  complex.remove_popable_blockers();
  std::cout << "complex complex" << complex.to_string();
  std::cout << std::endl << std::endl;

  BOOST_CHECK(complex.num_blockers() == 1);

  // test 2
  complex.clear();
  build_complete(4, complex);
  complex.add_vertex();
  complex.add_vertex();
  complex.add_edge_without_blockers(Vertex_handle(3), Vertex_handle(4));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(4));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(5));
  complex.add_edge_without_blockers(Vertex_handle(3), Vertex_handle(5));
  complex.add_edge_without_blockers(Vertex_handle(4), Vertex_handle(5));
  sigma1 = Simplex(Vertex_handle(1), Vertex_handle(2), Vertex_handle(3));
  sigma2 = Simplex(Vertex_handle(2), Vertex_handle(3), Vertex_handle(4));

  complex.add_blocker(sigma1);
  complex.add_blocker(sigma2);
  std::cout << "complex complex" << complex.to_string();
  std::cout << std::endl << std::endl;
  std::cout << "complex.RemovePopableBlockers();" << std::endl;
  complex.remove_popable_blockers();
  std::cout << "complex complex" << complex.to_string();

  std::cout << std::endl << std::endl;
  BOOST_CHECK(complex.num_blockers() == 0);
}
