/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Toplex_map.h>

#include <iostream>
#include <utility>  // for pair
#include <vector>
#include <cassert>

int main(int argc, char* const argv[]) {
  using Simplex = Gudhi::Toplex_map::Simplex;
  Simplex sigma1 = {1, 2, 3};
  Simplex sigma2 = {2, 3, 4, 5};

  Gudhi::Toplex_map tm;
  tm.insert_simplex(sigma1);
  tm.insert_simplex(sigma2);

  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*   /X\5/       */
  /*  o---o        */
  /*  1   3        */

  std::cout << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal cofaces
  Simplex sigma3 = {2, 3};
  std::cout << "Maximal cofaces of {2, 3} are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_cofaces(sigma3, 2)) {
    for (auto v : *simplex_ptr) {
      std::cout << v << ", ";
    }
    std::cout << std::endl;
  }

  // Browse maximal simplices
  std::cout << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::cout << v << ", ";
    }
    std::cout << std::endl;
  }

  Simplex sigma4 = {1, 3};
  assert(tm.membership(sigma4));

  Gudhi::Toplex_map::Vertex v = tm.contraction(1, 3);
  std::cout << "After contraction(1, 3) - " << v << std::endl;
  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*     \5/       */
  /*      o        */
  /*      3        */
  std::cout << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal simplices
  std::cout << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::cout << v << ", ";
    }
    std::cout << std::endl;
  }

  Simplex sigma5 = {3, 4};
  assert(tm.membership(sigma5));

  v = tm.contraction(3, 4);
  std::cout << "After contraction(3, 4) - " << v << std::endl;
  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*     \X/       */
  /*      o        */
  /*      5        */
  std::cout << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal simplices
  std::cout << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::cout << v << ", ";
    }
    std::cout << std::endl;
  }

  tm.insert_simplex(sigma1);
  tm.insert_simplex(sigma2);
  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*   /X\5/       */
  /*  o---o        */
  /*  1   3        */
  tm.remove_simplex(sigma1);

  std::cout << "After remove_simplex(1, 2, 3)" << std::endl;
  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*   / \5/       */
  /*  o---o        */
  /*  1   3        */
  std::cout << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal simplices
  std::cout << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::cout << v << ", ";
    }
    std::cout << std::endl;
  }

  tm.remove_vertex(1);

  std::cout << "After remove_vertex(1)" << std::endl;
  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*     \5/       */
  /*      o        */
  /*      3        */
  std::cout << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal simplices
  std::cout << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::cout << v << ", ";
    }
    std::cout << std::endl;
  }

  return 0;
}
