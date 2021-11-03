/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "permutahedral_representation"
#include <boost/test/unit_test.hpp>

#include <gudhi/Permutahedral_representation.h>

BOOST_AUTO_TEST_CASE(permutahedral_representation) {
  typedef std::vector<int> Vertex;
  typedef std::vector<std::size_t> Part;
  typedef std::vector<Part> Partition;
  typedef Gudhi::coxeter_triangulation::Permutahedral_representation<Vertex, Partition> Simplex_handle;
  Vertex v0(10, 0);
  Partition omega = {Part({5}), Part({2}), Part({3, 7}), Part({4, 9}), Part({0, 6, 8}), Part({1, 10})};
  Simplex_handle s(v0, omega);

  // Dimension check
  BOOST_CHECK(s.dimension() == 5);

  // Vertex number check
  std::vector<Vertex> vertices;
  for (auto& v : s.vertex_range()) vertices.push_back(v);
  BOOST_CHECK(vertices.size() == 6);

  // Facet number check
  std::vector<Simplex_handle> facets;
  for (auto& f : s.facet_range()) facets.push_back(f);
  BOOST_CHECK(facets.size() == 6);

  // Face of dim 3 number check
  std::vector<Simplex_handle> faces3;
  for (auto& f : s.face_range(3)) faces3.push_back(f);
  BOOST_CHECK(faces3.size() == 15);

  // Cofacet number check
  std::vector<Simplex_handle> cofacets;
  for (auto& f : s.cofacet_range()) cofacets.push_back(f);
  BOOST_CHECK(cofacets.size() == 12);

  // Is face check
  Vertex v1(10, 0);
  Partition omega1 = {Part({5}), Part({0, 1, 2, 3, 4, 6, 7, 8, 9, 10})};
  Simplex_handle s1(v1, omega1);
  Vertex v2(10, 0);
  v2[1] = -1;
  Partition omega2 = {Part({1}), Part({5}), Part({2}), Part({3, 7}), Part({4, 9}), Part({0, 6, 8}), Part({10})};
  Simplex_handle s2(v2, omega2);
  BOOST_CHECK(s.is_face_of(s));
  BOOST_CHECK(s1.is_face_of(s));
  BOOST_CHECK(!s2.is_face_of(s));
  BOOST_CHECK(s.is_face_of(s2));
}
