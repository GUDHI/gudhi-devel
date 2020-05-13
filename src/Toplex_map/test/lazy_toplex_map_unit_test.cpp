/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author:       François Godi, Vincent Rouvreau
 *
 *    Copyright (C) 2018  INRIA
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <vector>
#include <gudhi/Lazy_toplex_map.h>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "lazy toplex map"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(toplex_map) {
  using Vertex = Gudhi::Lazy_toplex_map::Vertex;

  Gudhi::Lazy_toplex_map tm;
  std::clog << "insert_simplex {1, 2, 3, 4}" << std::endl;
  std::vector<Vertex> sigma1 = {1, 2, 3, 4};
  tm.insert_simplex(sigma1);
  std::clog << "insert_simplex {5, 2, 3, 6}" << std::endl;
  std::vector<Vertex> sigma2 = {5, 2, 3, 6};
  tm.insert_simplex(sigma2);
  std::clog << "insert_simplex {5}" << std::endl;
  std::vector<Vertex> sigma3 = {5};
  tm.insert_simplex(sigma3);
  std::clog << "insert_simplex {4, 5, 3}" << std::endl;
  std::vector<Vertex> sigma6 = {4, 5, 3};
  tm.insert_simplex(sigma6);
  std::clog << "insert_simplex {4, 5, 9}" << std::endl;
  std::vector<Vertex> sigma7 = {4, 5, 9};
  tm.insert_simplex(sigma7);

  std::clog << "num_maximal_simplices = " << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 5);

  std::vector<Vertex> sigma4 = {5, 2, 3};
  std::vector<Vertex> sigma5 = {5, 2, 7};
  BOOST_CHECK(tm.membership(sigma4));
  BOOST_CHECK(!tm.membership(sigma5));
  std::clog << "insert_simplex {5, 2, 7}" << std::endl;
  tm.insert_simplex(sigma5);

  std::clog << "num_maximal_simplices = " << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 6);

  BOOST_CHECK(tm.membership(sigma5));

  std::clog << "contraction(4,5)" << std::endl;
  auto r = tm.contraction(4, 5);
  std::clog << "r=" << r << std::endl;
  BOOST_CHECK(r == 5);

  std::clog << "num_maximal_simplices = " << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 6);

  std::vector<Vertex> sigma8 = {1, 2, 3};
  std::vector<Vertex> sigma9 = {2, 7};

  sigma8.emplace_back(r);
  sigma9.emplace_back(r);
  BOOST_CHECK(!tm.membership(sigma6));
  BOOST_CHECK(tm.membership(sigma8));
  BOOST_CHECK(tm.membership(sigma9));

  std::clog << "remove_simplex({2, 7, r = 5})" << std::endl;
  tm.remove_simplex(sigma9);
  BOOST_CHECK(!tm.membership(sigma9));

  std::clog << "num_maximal_simplices = " << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 8);

  // {2, 7, 5} is removed, but verify its edges are still there
  std::vector<Vertex> edge = {2, 7};
  BOOST_CHECK(tm.membership(edge));
  edge = {2, 5};
  BOOST_CHECK(tm.membership(edge));
  edge = {7, 5};
  BOOST_CHECK(tm.membership(edge));
}

BOOST_AUTO_TEST_CASE(toplex_map_empty_toplex) {
  using Vertex = Gudhi::Lazy_toplex_map::Vertex;

  Gudhi::Lazy_toplex_map tm;
  std::clog << "num_maximal_simplices = " << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 0);
  std::clog << "num_vertices = " << tm.num_vertices() << std::endl;
  BOOST_CHECK(tm.num_vertices() == 0);

  std::clog << "Check an empty simplex is a member." << std::endl;
  std::vector<Vertex> empty_sigma = {};
  BOOST_CHECK(tm.membership(empty_sigma));

  std::clog << "Check the edge 2,7 is not a member." << std::endl;
  std::vector<Vertex> edge = {2, 7};
  BOOST_CHECK(!tm.membership(edge));

  std::clog << "Insert an empty simplex." << std::endl;
  tm.insert_simplex(empty_sigma);

  std::clog << "num_maximal_simplices = " << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 0);
  std::clog << "num_vertices = " << tm.num_vertices() << std::endl;
  BOOST_CHECK(tm.num_vertices() == 0);

  std::clog << "Check an empty simplex is a member." << std::endl;
  BOOST_CHECK(tm.membership(empty_sigma));
  std::clog << "Check the edge 2,7 is not a member." << std::endl;
  BOOST_CHECK(!tm.membership(edge));

  std::clog << "Insert edge 2,7." << std::endl;
  tm.insert_simplex(edge);

  std::clog << "num_maximal_simplices = " << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 1);
  std::clog << "num_vertices = " << tm.num_vertices() << std::endl;
  BOOST_CHECK(tm.num_vertices() == 2);

  std::clog << "Check an empty simplex is a member." << std::endl;
  BOOST_CHECK(tm.membership(empty_sigma));
  std::clog << "Check the edge 2,7 is a member." << std::endl;
  BOOST_CHECK(tm.membership(edge));

  std::clog << "contraction(2,7)" << std::endl;
  auto r = tm.contraction(2, 7);
  std::clog << "r=" << r << std::endl;
  BOOST_CHECK(r == 7);

  std::clog << "num_maximal_simplices = " << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 1);
  std::clog << "num_vertices = " << tm.num_vertices() << std::endl;
  BOOST_CHECK(tm.num_vertices() == 1);

  std::clog << "Check an empty simplex is a member." << std::endl;
  BOOST_CHECK(tm.membership(empty_sigma));
  std::clog << "Check the edge 2,7 is not a member." << std::endl;
  BOOST_CHECK(!tm.membership(edge));

  std::clog << "Remove the vertex 7." << std::endl;
  std::vector<Vertex> vertex = {7};
  tm.remove_simplex(vertex);

  std::clog << "num_maximal_simplices = " << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 0);
  std::clog << "num_vertices = " << tm.num_vertices() << std::endl;
  BOOST_CHECK(tm.num_vertices() == 0);

  std::clog << "Check an empty simplex is a member." << std::endl;
  BOOST_CHECK(tm.membership(empty_sigma));
  std::clog << "Check the edge 2,7 is not a member." << std::endl;
  BOOST_CHECK(!tm.membership(edge));
}
