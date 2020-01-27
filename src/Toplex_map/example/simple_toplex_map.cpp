/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
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

  std::clog << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal cofaces
  Simplex sigma3 = {2, 3};
  std::clog << "Maximal cofaces of {2, 3} are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_cofaces(sigma3, 2)) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
  }

  // Browse maximal simplices
  std::clog << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
  }

  Simplex sigma4 = {1, 3};
  assert(tm.membership(sigma4));

  Gudhi::Toplex_map::Vertex v = tm.contraction(1, 3);
  std::clog << "After contraction(1, 3) - " << v << std::endl;
  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*     \5/       */
  /*      o        */
  /*      3        */
  std::clog << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal simplices
  std::clog << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
  }

  Simplex sigma5 = {3, 4};
  assert(tm.membership(sigma5));

  v = tm.contraction(3, 4);
  std::clog << "After contraction(3, 4) - " << v << std::endl;
  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*     \X/       */
  /*      o        */
  /*      5        */
  std::clog << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal simplices
  std::clog << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
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

  std::clog << "After remove_simplex(1, 2, 3)" << std::endl;
  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*   / \5/       */
  /*  o---o        */
  /*  1   3        */
  std::clog << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal simplices
  std::clog << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
  }

  tm.remove_vertex(1);

  std::clog << "After remove_vertex(1)" << std::endl;
  /* Simplex is:   */
  /*    2   4      */
  /*    o---o      */
  /*     \5/       */
  /*      o        */
  /*      3        */
  std::clog << "num max simplices = " << tm.num_maximal_simplices() << " - num vertices = " << tm.num_vertices()
            << std::endl;

  // Browse maximal simplices
  std::clog << "Maximal simplices are :" << std::endl;
  for (auto simplex_ptr : tm.maximal_simplices()) {
    for (auto v : *simplex_ptr) {
      std::clog << v << ", ";
    }
    std::clog << std::endl;
  }

  return 0;
}
