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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "skeleton_blocker_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Skeleton_blocker.h>

template<typename ComplexType> class Skeleton_blocker_sub_complex;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_simple_traits Traits;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_complex<Traits> Complex;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_link_complex<Complex> Skeleton_blocker_link_complex;
typedef Complex::Vertex_handle Vertex_handle;
typedef Complex::Root_vertex_handle Root_vertex_handle;
typedef Complex::Simplex Simplex;
typedef Complex::Root_simplex_handle Root_simplex_handle;
typedef Simplex::Simplex_vertex_const_iterator Simplex_vertex_const_iterator;
typedef Complex::Edge_handle Edge_handle;

bool assert_vertex(Complex &complex, Vertex_handle v) {
  //assert(complex.contains(v));
  return complex.contains(static_cast<Simplex> (v));
}

bool assert_simplex(Complex &complex, Root_vertex_handle a, Root_vertex_handle b, Root_vertex_handle c) {
  return true;
  //	AddressSimplex simplex((a),(b),(c));
  //	return complex.contains(&simplex);
}

// true iff the blocker (a,b,c) is in complex

bool assert_blocker(Complex &complex, Root_vertex_handle a, Root_vertex_handle b, Root_vertex_handle c) {
  return true;
  //return complex.contains_blocker((a),(b),(c));
}

// true iff the blocker (a,b,c,d) is in complex

bool assert_blocker(Complex &complex, Root_vertex_handle a, Root_vertex_handle b, Root_vertex_handle c, Root_vertex_handle d) {
  return true;
  //Simplex blocker (a,b,c,d);
  //return complex.contains_blocker(&blocker);
}

void build_complete(int n, Complex& complex) {
  complex.clear();
  for (int i = 0; i < n; i++)
    complex.add_vertex();

  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_simplex) {
  Simplex simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2), Vertex_handle(3));
  BOOST_CHECK(simplex.dimension() == 3);
}

BOOST_AUTO_TEST_CASE(test_skeleton_num_simplices) {
  int n = 4;
  Complex complex;
  build_complete(n, complex);
  size_t sum = 0;
  for (int i = 0; i < n; i++) {
    sum += complex.num_simplices(i);
  }
  BOOST_CHECK(complex.num_vertices() == n);
  BOOST_CHECK(complex.num_edges() == 6);
  BOOST_CHECK(sum == 15);
  BOOST_CHECK(complex.num_simplices() == 15);
}


BOOST_AUTO_TEST_CASE(test_skeleton_iterator_vertices1) {
  int n = 10;
  Complex complex(10);
  std::cout << "complex.num_vertices():" << complex.num_vertices() << std::endl;
  int num_vertex_seen = 0;
  for (auto vi : complex.vertex_range()) {
    std::cout << "vertex:" << vi << std::endl;
    ++num_vertex_seen;
  }
  BOOST_CHECK(num_vertex_seen == n);
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_vertices2) {
  int n = 10;
  Complex complex;
  build_complete(10, complex);
  std::cout << "complex.num_vertices():" << complex.num_vertices() << std::endl;
  std::cout << "complex.num_edges():" << complex.num_edges() << std::endl;
  int num_vertex_seen = 0;
  for (auto vi : complex.vertex_range(Vertex_handle(2))) {
    std::cout << "vertex:" << vi << std::endl;
    ++num_vertex_seen;
  }
  std::cout << "num_vertex_seen:" << num_vertex_seen << std::endl;
  BOOST_CHECK(num_vertex_seen == (n -1));
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_edge) {
  const int n = 10;
  Complex complex(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
  complex.remove_edge(Vertex_handle(2), Vertex_handle(3));
  complex.remove_edge(Vertex_handle(3), Vertex_handle(5));
  std::cout << "complex.num_edges():" << complex.num_edges() << std::endl;
  int num_edges_seen = 0;
  for (auto edge : complex.edge_range()) {
    std::cout << "edge :" << complex[edge] << std::endl;
    ++num_edges_seen;
  }

  BOOST_CHECK(num_edges_seen == n * (n - 1) / 2 - 2);
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_edge2) {
  const int n = 10;
  Complex complex(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
  complex.remove_edge(Vertex_handle(2), Vertex_handle(3));
  complex.remove_edge(Vertex_handle(3), Vertex_handle(5));
  std::cout << "complex.num_edges():" << complex.num_edges() << std::endl;
  int num_neigbors_seen = 0;
  for (auto neighbor : complex.vertex_range(Vertex_handle(2))) {
    std::cout << "neighbor" << neighbor << std::endl;
    ++num_neigbors_seen;
  }
  BOOST_CHECK(num_neigbors_seen == 8);
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_triangles) {
  const int n = 7;
  Complex complex(n);
  //create a "ring" around '0'
  for (int i = 1; i < n; i++)
    complex.add_edge_without_blockers(Vertex_handle(0), Vertex_handle(i));
  for (int i = 1; i < n - 1; i++)
    complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(i + 1));
  complex.add_edge_without_blockers(Vertex_handle(1), Vertex_handle(6));

  std::cout << complex.to_string() << std::endl;

  int num_triangles_seen = 0;
  //for (auto t : complex.triangle_range(5)){
  for (auto t : complex.triangle_range(Vertex_handle(5))) {
    ++num_triangles_seen;
  }
  BOOST_CHECK(num_triangles_seen == 2);

  num_triangles_seen = 0;
  for (auto t : complex.triangle_range(Vertex_handle(0))) {
    ++num_triangles_seen;
  }
  BOOST_CHECK(num_triangles_seen == 6);

  // we now add another triangle
  complex.add_vertex();
  complex.add_edge_without_blockers(Vertex_handle(4), Vertex_handle(7));
  complex.add_edge_without_blockers(Vertex_handle(3), Vertex_handle(7));
  complex.add_blocker(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(6)));
  num_triangles_seen = 0;

  num_triangles_seen = 0;
  for (auto t : complex.triangle_range()) {
    ++num_triangles_seen;
  }
  BOOST_CHECK(num_triangles_seen == 6);
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_simplices) {
  Complex complex(6);
  complex.add_edge_without_blockers(Vertex_handle(0), Vertex_handle(1));
  complex.add_edge_without_blockers(Vertex_handle(1), Vertex_handle(2));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(0));
  complex.add_edge_without_blockers(Vertex_handle(1), Vertex_handle(3));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(3));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(5));
  complex.add_edge_without_blockers(Vertex_handle(3), Vertex_handle(5));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(4));
  complex.add_edge_without_blockers(Vertex_handle(4), Vertex_handle(5));
  complex.add_edge_without_blockers(Vertex_handle(3), Vertex_handle(4));

  complex.add_blocker(Simplex(Vertex_handle(2), Vertex_handle(3), Vertex_handle(4), Vertex_handle(5)));

  std::map<Vertex_handle, unsigned> expected_num_simplices;

  expected_num_simplices[Vertex_handle(0)] = 4;
  expected_num_simplices[Vertex_handle(1)] = 6;
  expected_num_simplices[Vertex_handle(2)] = 11;
  expected_num_simplices[Vertex_handle(3)] = 9;
  expected_num_simplices[Vertex_handle(4)] = 7;
  expected_num_simplices[Vertex_handle(5)] = 7;

  for (auto pair : expected_num_simplices) {
    std::cout << "found list: ";
    unsigned num_simplices_around = 0;
    for (const auto& simplex : complex.star_simplex_range(pair.first)) {
      simplex.dimension();
      std::cout << simplex << " - ";
      ++num_simplices_around;
    }

    BOOST_CHECK(num_simplices_around == pair.second);

    std::cout << std::endl << "current vertex:" << pair.first << " - ";
    std::cout << "expected_num_simplices:" << pair.second << " - ";
    std::cout << "found:" << num_simplices_around << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_simplices2) {
  Complex complex(2);
  complex.add_edge_without_blockers(Vertex_handle(0), Vertex_handle(1));

  // Check there is no triangle
  BOOST_CHECK(std::distance(complex.triangle_range().begin(), complex.triangle_range().end()) == 0);

  // Star(0) is [{0},{0,1}]
  BOOST_CHECK(std::distance(complex.star_simplex_range(Vertex_handle(0)).begin(),
                            complex.star_simplex_range(Vertex_handle(0)).end()) == 2);
  
  // No blocker
  BOOST_CHECK(std::distance(complex.blocker_range(Vertex_handle(0)).begin(),
                            complex.blocker_range(Vertex_handle(0)).end()) == 0);
  
  // Complex is [{0},{0,1},{1}]
  BOOST_CHECK(std::distance(complex.complex_simplex_range().begin(),
                            complex.complex_simplex_range().end()) == 3);
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_simplices3) {
  Complex complex(3);
  complex.add_edge_without_blockers(Vertex_handle(0), Vertex_handle(1));
  complex.add_edge_without_blockers(Vertex_handle(1), Vertex_handle(2));
  complex.add_edge_without_blockers(Vertex_handle(2), Vertex_handle(0));
  complex.add_blocker(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2)));

  // Check there is no triangle
  BOOST_CHECK(std::distance(complex.triangle_range().begin(), complex.triangle_range().end()) == 0);

  // Star(0) is [{0},{0,1},{0,2}]
  BOOST_CHECK(std::distance(complex.star_simplex_range(Vertex_handle(0)).begin(),
                            complex.star_simplex_range(Vertex_handle(0)).end()) == 3);

  // blocker(0) is [{0,1,2}]
  BOOST_CHECK(std::distance(complex.blocker_range(Vertex_handle(0)).begin(),
                            complex.blocker_range(Vertex_handle(0)).end()) == 1);
  
  // Complex is [{0},{0,1},{0,2},{1},{1,2},{2}]
  BOOST_CHECK(std::distance(complex.complex_simplex_range().begin(),
                            complex.complex_simplex_range().end()) == 6);
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_simplices4) {
  Complex empty_complex;
  for (auto v : empty_complex.vertex_range()) {
    std::cout << v;
    BOOST_CHECK(false);
  }
  for (auto e : empty_complex.edge_range()) {
    std::cout << e;
    BOOST_CHECK(false);
  }
  for (auto t : empty_complex.triangle_range()) {
    std::cout << t;
    BOOST_CHECK(false);
  }
  for (auto s : empty_complex.complex_simplex_range()) {
    std::cout << s;
    BOOST_CHECK(false);
  }
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_coboundary) {
  Complex c;
  build_complete(4, c);
  c.remove_edge(Vertex_handle(1), Vertex_handle(3));
  std::cout << c.to_string();
  Simplex s02(Vertex_handle(0), Vertex_handle(2));
  int n = 0;
  std::set<Simplex> expected;
  expected.insert(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2)));
  expected.insert(Simplex(Vertex_handle(0), Vertex_handle(2), Vertex_handle(3)));
  for (const auto & s : c.coboundary_range(s02)) {
    BOOST_CHECK(expected.find(s) != expected.end());
    ++n;
  }
  BOOST_CHECK(n == 2);
}

template<typename Map>
auto blocker_range(Map map) -> decltype(map | boost::adaptors::map_values) {
  return map | boost::adaptors::map_values;
}

BOOST_AUTO_TEST_CASE(test_skeleton_iterator_blockers) {
  Complex complex;
  Simplex alpha;
  Simplex vertex_set_expected;
  // Build the complexes
  for (int i = 0; i < 20; i++) {
    complex.add_vertex();
  }
  for (int i = 10; i < 15; i++) {
    for (int j = i + 1; j < 15; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
  }

  std::vector<Simplex> myBlockers;
  myBlockers.push_back(Simplex(Vertex_handle(10), Vertex_handle(11), Vertex_handle(12)));
  myBlockers.push_back(Simplex(Vertex_handle(2), Vertex_handle(1), Vertex_handle(10)));
  myBlockers.push_back(Simplex(Vertex_handle(10), Vertex_handle(9), Vertex_handle(15)));
  myBlockers.push_back(Simplex(Vertex_handle(1), Vertex_handle(9), Vertex_handle(8)));

  for (auto blocker : myBlockers)
    complex.add_blocker(blocker);

  int num_blockers = 0;
  for (auto blockers : complex.blocker_range(Vertex_handle(10))) {
    // Only the first 3 blockers contain vertex 10
    BOOST_CHECK(*blockers == myBlockers[num_blockers]);
    num_blockers++;
  }
  BOOST_CHECK(num_blockers == 3);

  num_blockers = 0;
  for (auto blockers : complex.blocker_range()) {
// If not windows - _WIN32 is for windows 32 and 64 bits
#ifndef _WIN32
    for (auto block_ptr = myBlockers.begin(); block_ptr < myBlockers.end(); block_ptr++)
      if (*block_ptr == *blockers)
        myBlockers.erase(block_ptr);
#endif
    num_blockers++;
  }
  BOOST_CHECK(num_blockers == 4);
// If not windows - _WIN32 is for windows 32 and 64 bits
#ifndef _WIN32
  BOOST_CHECK(myBlockers.empty());
#endif
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_link0) {
  enum { a, b, c, d, n };
  Complex complex(n);
  complex.add_edge_without_blockers(Vertex_handle(b), Vertex_handle(c));
  complex.add_edge_without_blockers(Vertex_handle(c), Vertex_handle(d));
  Simplex alpha = Simplex(Vertex_handle(c));
  Skeleton_blocker_link_complex L(complex, alpha);

  auto L2 = complex.link(alpha);
  BOOST_CHECK(L == L2);

  std::cout << L.to_string();

  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(b))));
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(d))));
  BOOST_CHECK(L.num_edges() == 0);
  BOOST_CHECK(L.num_blockers() == 0);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_link1) {
  Complex complex;

  // Build the complexes
  for (int i = 0; i < 20; i++) {
    complex.add_vertex();
  }
  for (int i = 10; i < 15; i++) {
    for (int j = i + 1; j < 15; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
  }
  Simplex alpha(Vertex_handle(12), Vertex_handle(14));
  Skeleton_blocker_link_complex L(complex, alpha);
  // Complexes built

  auto L2 = complex.link(alpha);
  BOOST_CHECK(L == L2);

  // verification
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(10))));
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(11))));
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(13))));
  BOOST_CHECK(L.num_edges() == 3);
  BOOST_CHECK(L.num_blockers() == 0);
  Root_simplex_handle simplex;
  simplex.add_vertex(Root_vertex_handle(10));
  simplex.add_vertex(Root_vertex_handle(11));
  simplex.add_vertex(Root_vertex_handle(13));
  BOOST_CHECK(L.get_simplex_address(simplex));
  BOOST_CHECK(L.contains(*(L.get_simplex_address(simplex))));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_link2) {
  Complex complex;

  Simplex alpha;
  Simplex vertex_set_expected;
  // Build the complexes
  for (int i = 0; i < 20; i++) {
    complex.add_vertex();
  }
  for (int i = 10; i < 15; i++) {
    for (int j = i + 1; j < 15; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
  }
  complex.add_blocker(Simplex(Vertex_handle(10), Vertex_handle(11), Vertex_handle(13)));
  alpha = Simplex(Vertex_handle(12), Vertex_handle(14));
  Skeleton_blocker_link_complex L(complex, alpha);
  // Complexes built

  // Print result
  std::cout << "complex complex" << complex.to_string();
  std::cout << std::endl << std::endl;
  std::cout << "L= Link_complex(" << alpha << ") : \n" << L.to_string();

  auto L2 = complex.link(alpha);
  BOOST_CHECK(L == L2);


  // verification
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(10))));
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(11))));
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(13))));
  BOOST_CHECK(L.num_edges() == 3);
  BOOST_CHECK(L.num_blockers() == 1);
  Root_simplex_handle simplex;
  simplex.add_vertex(Root_vertex_handle(10));
  simplex.add_vertex(Root_vertex_handle(11));
  simplex.add_vertex(Root_vertex_handle(13));
  BOOST_CHECK(L.contains_blocker(*(L.get_simplex_address(simplex))));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_link3) {
  Complex complex;

  Simplex alpha;
  Simplex vertex_set_expected;
  // Build the complexes
  for (int i = 0; i < 20; i++) {
    complex.add_vertex();
  }
  for (int i = 10; i < 15; i++) {
    for (int j = i + 1; j < 15; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
  }
  complex.add_blocker(Simplex(Vertex_handle(10), Vertex_handle(11), Vertex_handle(12)));
  alpha = Simplex(Vertex_handle(12), Vertex_handle(14));
  Skeleton_blocker_link_complex L(complex, alpha);
  // Complexes built

  // Print result
  std::cout << "complex complex" << complex.to_string();
  std::cout << std::endl << std::endl;
  std::cout << "L= Link_complex(" << alpha << ") : \n" << L.to_string();

  auto L2 = complex.link(alpha);
  BOOST_CHECK(L == L2);

  // verification
  BOOST_CHECK(L.contains(static_cast<Simplex> (*L.get_address(Root_vertex_handle(10)))));
  BOOST_CHECK(L.contains(static_cast<Simplex> (*L.get_address(Root_vertex_handle(11)))));
  BOOST_CHECK(L.contains(static_cast<Simplex> (*L.get_address(Root_vertex_handle(13)))));
  BOOST_CHECK(L.num_edges() == 2);
  BOOST_CHECK(L.contains_edge(*L.get_address(Root_vertex_handle(10)), *L.get_address(Root_vertex_handle(13))));
  BOOST_CHECK(L.contains_edge(*L.get_address(Root_vertex_handle(13)), *L.get_address(Root_vertex_handle(11))));
  BOOST_CHECK(L.num_blockers() == 0);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_link4) {
  Complex complex;

  // Build the complexes
  for (int i = 0; i < 20; i++) {
    complex.add_vertex();
  }
  for (int i = 10; i < 15; i++) {
    for (int j = i + 1; j < 15; j++)
      complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(j));
  }
  complex.add_blocker(Simplex(Vertex_handle(10), Vertex_handle(11), Vertex_handle(12), Vertex_handle(13)));
  Simplex alpha(Vertex_handle(12), Vertex_handle(14));
  Skeleton_blocker_link_complex L(complex, alpha);
  // Complexes built

  // verification
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(10))));
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(11))));
  BOOST_CHECK(L.contains_vertex(*L.get_address(Root_vertex_handle(13))));
  BOOST_CHECK(L.num_edges() == 3);
  BOOST_CHECK(L.num_blockers() == 1);
  Root_simplex_handle simplex;
  simplex.add_vertex(Root_vertex_handle(10));
  simplex.add_vertex(Root_vertex_handle(11));
  simplex.add_vertex(Root_vertex_handle(13));
  BOOST_CHECK(L.contains_blocker(*(L.get_simplex_address(simplex))));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_link5) {
  Complex complex(0);
  // Build the complexes
  build_complete(4, complex);
  complex.add_blocker(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2), Vertex_handle(3)));

  Simplex alpha(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2));
  Skeleton_blocker_link_complex L(complex, alpha);
  // Complexes built

  // Print result
  std::cout << "Complex: " << complex.to_string()<< std::endl << std::endl;
  std::cout << "Link: " << L.to_string() << std::endl;

  // verification
  BOOST_CHECK(L.num_vertices() == 0);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_link6) {
  Complex complex(0);
  // Build the complexes
  build_complete(4, complex);
  complex.add_blocker(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2)));

  Simplex alpha(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2));

  Skeleton_blocker_link_complex link_blocker_alpha;

  build_link_of_blocker(complex, alpha, link_blocker_alpha);

  // Print result
  std::cout << "Complex: " << complex.to_string()<< std::endl << std::endl;
  std::cout << "Link: " << link_blocker_alpha.to_string() << std::endl;

  // verification
  BOOST_CHECK(link_blocker_alpha.num_vertices() == 1);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_link7) {
  Complex complex(0);
  // Build the complexes
  build_complete(6, complex);
  complex.add_vertex();
  complex.add_vertex();
  for (int i = 3; i < 6; ++i) {
    complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(6));
    complex.add_edge_without_blockers(Vertex_handle(i), Vertex_handle(7));
  }
  complex.add_edge_without_blockers(Vertex_handle(6), Vertex_handle(7));
  complex.add_blocker(Simplex(Vertex_handle(0), Vertex_handle(1), Vertex_handle(2)));
  complex.add_blocker(Simplex(Vertex_handle(3), Vertex_handle(4), Vertex_handle(5)));

  Simplex alpha(Vertex_handle(3), Vertex_handle(4), Vertex_handle(5));

  Skeleton_blocker_link_complex link_blocker_alpha;

  build_link_of_blocker(complex, alpha, link_blocker_alpha);

  //the result should be the edge {6,7} plus the blocker {0,1,2}

  // Print result
  std::cout << "Complex: " << complex.to_string()<< std::endl << std::endl;
  std::cout << "Link: " << link_blocker_alpha.to_string() << std::endl;

  Skeleton_blocker_link_complex link_blocker_alpha_cpy = link_blocker_alpha;

  std::cout << "Link copy: " << link_blocker_alpha_cpy.to_string() << std::endl;

  BOOST_CHECK(link_blocker_alpha.num_vertices() == link_blocker_alpha_cpy.num_vertices());
  BOOST_CHECK(link_blocker_alpha.num_blockers() == link_blocker_alpha_cpy.num_blockers());
  BOOST_CHECK(link_blocker_alpha.num_edges() == link_blocker_alpha_cpy.num_edges());
  BOOST_CHECK((link_blocker_alpha.num_blockers() == link_blocker_alpha_cpy.num_blockers()));

  // verification
  BOOST_CHECK(link_blocker_alpha.num_vertices() == 5);
  BOOST_CHECK(link_blocker_alpha.num_edges() == 4);
  BOOST_CHECK(link_blocker_alpha.num_blockers() == 1);
}

template<typename SimplexHandle>
void add_triangle_edges(int a, int b, int c, std::list<SimplexHandle>& simplices) {
  typedef SimplexHandle Simplex;
  typedef typename SimplexHandle::Vertex_handle Vertex_handle;

  simplices.push_back(Simplex(Vertex_handle(a), Vertex_handle(b)));
  simplices.push_back(Simplex(Vertex_handle(b), Vertex_handle(c)));
  simplices.push_back(Simplex(Vertex_handle(c), Vertex_handle(a)));
}

template<typename SimplexHandle>
void add_triangle(int a, int b, int c, std::list<SimplexHandle>& simplices) {
  typedef SimplexHandle Simplex;
  typedef typename SimplexHandle::Vertex_handle Vertex_handle;
  simplices.push_back(Simplex(Vertex_handle(a), Vertex_handle(b), Vertex_handle(c)));
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_constructor) {
  std::list <Simplex> simplices;

  simplices.push_back(Simplex(Vertex_handle(0)));
  simplices.push_back(Simplex(Vertex_handle(1)));
  simplices.push_back(Simplex(Vertex_handle(2)));
  simplices.push_back(Simplex(Vertex_handle(3)));
  simplices.push_back(Simplex(Vertex_handle(4)));
  simplices.push_back(Simplex(Vertex_handle(5)));

  simplices.push_back(Simplex(Vertex_handle(3), Vertex_handle(5)));

  add_triangle_edges(0, 1, 5, simplices);
  add_triangle_edges(1, 2, 3, simplices);
  add_triangle_edges(2, 3, 4, simplices);
  add_triangle_edges(1, 3, 4, simplices);
  add_triangle_edges(1, 2, 4, simplices);

  add_triangle(0, 1, 5, simplices);
  add_triangle(1, 2, 3, simplices);
  add_triangle(1, 3, 4, simplices);
  add_triangle(1, 2, 4, simplices);
  add_triangle(2, 3, 4, simplices);

  Complex complex(simplices.begin(), simplices.end());

  std::cout << "Constructor 1:\n" << complex.to_string();

  BOOST_CHECK(complex.num_vertices() == 6);
  BOOST_CHECK(complex.num_edges() == 10);
  BOOST_CHECK(complex.num_blockers() == 2);
}

std::list<Simplex> subfaces(Simplex top_face) {
  std::list<Simplex> res;
  if (top_face.dimension() == -1) return res;
  if (top_face.dimension() == 0) {
    res.push_back(top_face);
    return res;
  } else {
    Vertex_handle first_vertex = top_face.first_vertex();
    top_face.remove_vertex(first_vertex);
    res = subfaces(top_face);
    std::list<Simplex> copy = res;
    for (auto& simplex : copy) {
      simplex.add_vertex(first_vertex);
    }
    res.push_back(Simplex(first_vertex));
    res.splice(res.end(), copy);
    return res;
  }
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_constructor2) {
  Simplex simplex;
  for (int i = 0; i < 5; ++i)
    simplex.add_vertex(static_cast<Vertex_handle> (i));

  std::list <Simplex> simplices(subfaces(simplex));
  simplices.remove(simplex);

  Complex complex(simplices.begin(), simplices.end());

  std::cout << "Constructor 2:\n" << complex.to_string();

  for (auto b : complex.const_blocker_range()) {
    std::cout << "b:" << b << std::endl;
  }

  BOOST_CHECK(complex.num_vertices() == 5);
  BOOST_CHECK(complex.num_edges() == 10);
  BOOST_CHECK(complex.num_blockers() == 1);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_constructor3) {
  typedef Vertex_handle Vh;
  typedef Simplex Sh;
  std::vector<Simplex> simplices;
  auto subf(subfaces(Sh(Vh(0), Vh(1), Vh(2))));
  subf.pop_back(); //remove max face -> now a blocker 012
  simplices.insert(simplices.begin(), subf.begin(), subf.end());

  Complex complex(simplices.begin(), simplices.end());

  std::cout << "Constructor 3:\n" << complex.to_string();

  BOOST_CHECK(complex.num_blockers() == 1);
  Sh expected_blocker(Vh(0), Vh(1), Vh(2));
  for (auto b : complex.const_blocker_range())
    BOOST_CHECK(*b == expected_blocker);


  BOOST_CHECK(complex.num_vertices() == 3);
  BOOST_CHECK(complex.num_blockers() == 1);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_constructor4) {
  typedef Vertex_handle Vh;
  typedef Simplex Sh;
  std::vector<Simplex> simplices;
  auto subf(subfaces(Sh(Vh(0), Vh(1), Vh(2), Vh(3))));
  simplices.insert(simplices.begin(), subf.begin(), subf.end());

  simplices.push_back(Sh(Vh(4)));
  simplices.push_back(Sh(Vh(4), Vh(1)));
  simplices.push_back(Sh(Vh(4), Vh(0)));

  Complex complex(simplices.begin(), simplices.end());

  std::cout << "Constructor 4:\n" << complex.to_string();
  BOOST_CHECK(complex.num_blockers() == 1);
  Sh expected_blocker(Vh(0), Vh(1), Vh(4));
  for (auto b : complex.const_blocker_range())
    BOOST_CHECK(*b == expected_blocker);

  BOOST_CHECK(complex.num_vertices() == 5);
  BOOST_CHECK(complex.num_blockers() == 1);
  BOOST_CHECK(complex.num_edges() == 8);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_constructor5) {
  typedef Vertex_handle Vh;
  typedef Simplex Sh;
  std::vector<Simplex> simplices;
  auto subf(subfaces(Sh(Vh(0), Vh(1), Vh(2))));
  simplices.insert(simplices.begin(), subf.begin(), subf.end());

  simplices.push_back(Sh(Vh(3)));
  simplices.push_back(Sh(Vh(3), Vh(1)));
  simplices.push_back(Sh(Vh(3), Vh(2)));
  simplices.push_back(Sh(Vh(4)));
  simplices.push_back(Sh(Vh(4), Vh(1)));
  simplices.push_back(Sh(Vh(4), Vh(0)));
  simplices.push_back(Sh(Vh(5)));
  simplices.push_back(Sh(Vh(5), Vh(2)));
  simplices.push_back(Sh(Vh(5), Vh(0)));

  Complex complex(simplices.begin(), simplices.end());

  std::cout << "Constructor 5:\n" << complex.to_string();

  BOOST_CHECK(complex.num_vertices() == 6);
  BOOST_CHECK(complex.num_blockers() == 3);
  BOOST_CHECK(complex.num_edges() == 9);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_constructor6) {
  typedef Vertex_handle Vh;
  typedef Simplex Sh;
  std::vector<Simplex> simplices;
  auto subf(subfaces(Sh(Vh(0), Vh(1), Vh(2), Vh(3))));
  for (auto s : subf) {
    Sh s1(Vh(0), Vh(1), Vh(2), Vh(3));
    Sh s2(Vh(1), Vh(2), Vh(3));
    if (s != s1 && s != s2) simplices.push_back(s);
  }

  Complex complex(simplices.begin(), simplices.end());

  std::cout << "Constructor 6:\n" << complex.to_string();

  BOOST_CHECK(complex.num_vertices() == 4);
  BOOST_CHECK(complex.num_blockers() == 1);
  BOOST_CHECK(complex.num_edges() == 6);
  Sh expected_blocker(Vh(1), Vh(2), Vh(3));
  for (auto b : complex.const_blocker_range())
    BOOST_CHECK(*b == expected_blocker);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_constructor7) {
  typedef Vertex_handle Vh;
  typedef Simplex Sh;
  std::vector<Simplex> simplices;
  simplices.push_back(Sh(Vh(0), Vh(1), Vh(2)));
  simplices.push_back(Sh(Vh(1), Vh(2), Vh(3)));
  simplices.push_back(Sh(Vh(3), Vh(0), Vh(2)));
  simplices.push_back(Sh(Vh(3), Vh(0), Vh(1)));

  //get complex from top faces
  Complex complex(Gudhi::skeleton_blocker::make_complex_from_top_faces<Complex>(simplices.begin(), simplices.end()));

  std::cout << "Constructor 7:\n" << complex.to_string();

  BOOST_CHECK(complex.num_vertices() == 4);
  BOOST_CHECK(complex.num_blockers() == 1);
  BOOST_CHECK(complex.num_edges() == 6);
  Sh expected_blocker(Vh(0), Vh(1), Vh(2), Vh(3));
  for (auto b : complex.const_blocker_range())
    BOOST_CHECK(*b == expected_blocker);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_complex_constructor8) {
  typedef Vertex_handle Vh;
  typedef Simplex Sh;
  std::vector<Simplex> simplices;
  simplices.push_back(Sh(Vh(0), Vh(1)));
  simplices.push_back(Sh(Vh(2), Vh(1)));
  simplices.push_back(Sh(Vh(0), Vh(2)));
  simplices.push_back(Sh(Vh(3), Vh(1)));
  simplices.push_back(Sh(Vh(2), Vh(3)));

  //get complex from top faces
  Complex complex(Gudhi::skeleton_blocker::make_complex_from_top_faces<Complex>(simplices.begin(), simplices.end()));

  std::cout << "Constructor 8:\n" << complex.to_string();

  BOOST_CHECK(complex.num_vertices() == 4);
  BOOST_CHECK(complex.num_blockers() == 2);
  BOOST_CHECK(complex.num_edges() == 5);
}
