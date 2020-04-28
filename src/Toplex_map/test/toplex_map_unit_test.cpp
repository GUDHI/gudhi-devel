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
#include <gudhi/Toplex_map.h>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "toplex map"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(toplex_map) {
  using Vertex = Gudhi::Toplex_map::Vertex;

  Gudhi::Toplex_map tm;
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

  std::clog << "num_maximal_simplices" << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 4);
  // Browse maximal simplices
  std::clog << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
    BOOST_CHECK(tm.maximality(*simplex_ptr));
  }

  BOOST_CHECK(tm.maximality(sigma1));
  BOOST_CHECK(tm.maximality(sigma2));
  BOOST_CHECK(!tm.maximality(sigma3));
  BOOST_CHECK(tm.maximality(sigma6));
  BOOST_CHECK(tm.maximality(sigma7));

  std::vector<Vertex> sigma4 = {5, 2, 3};
  std::vector<Vertex> sigma5 = {5, 2, 7};
  BOOST_CHECK(tm.membership(sigma4));
  BOOST_CHECK(!tm.membership(sigma5));
  std::clog << "insert_simplex {5, 2, 7}" << std::endl;
  tm.insert_simplex(sigma5);

  std::clog << "num_maximal_simplices" << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 5);
  // Browse maximal simplices
  std::clog << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
    BOOST_CHECK(tm.maximality(*simplex_ptr));
  }

  BOOST_CHECK(tm.membership(sigma5));

  std::clog << "contraction(4,5)" << std::endl;
  auto r = tm.contraction(4, 5);
  std::clog << "r=" << r << std::endl;
  BOOST_CHECK(r == 5);

  std::clog << "num_maximal_simplices" << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 4);
  // Browse maximal simplices
  std::clog << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
    BOOST_CHECK(tm.maximality(*simplex_ptr));
  }

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

  std::clog << "num_maximal_simplices" << tm.num_maximal_simplices() << std::endl;
  BOOST_CHECK(tm.num_maximal_simplices() == 5);
  // Browse maximal simplices
  std::clog << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
    BOOST_CHECK(tm.maximality(*simplex_ptr));
  }
  // {2, 7, 5} is removed, but verify its edges are still there
  std::vector<Vertex> edge = {2, 7};
  BOOST_CHECK(tm.membership(edge));
  edge = {2, 5};
  BOOST_CHECK(tm.membership(edge));
  edge = {7, 5};
  BOOST_CHECK(tm.membership(edge));
}
