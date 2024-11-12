/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <string>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_constness"
#include <boost/test/unit_test.hpp>

#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

template<typename Simplex_tree>
void test_simplex_tree_constness(const Simplex_tree& const_stree) {
  Simplex_tree copy_operator = const_stree;
  Simplex_tree copy_ctor(const_stree);
  BOOST_CHECK(copy_operator == const_stree);
  BOOST_CHECK(const_stree == copy_ctor);

  std::clog << "********************************************************************\n";
  std::clog << "* The complex contains " << const_stree.num_simplices() << " simplices, " << const_stree.num_vertices();
  std::clog << " vertices - dimension " << const_stree.dimension() << "\n";

  std::clog << "* num_simplices_by_dimension = ";
  for (auto simplices_by_dim: const_stree.num_simplices_by_dimension())
    std::clog << simplices_by_dim << ",";
  std::clog << "\n";

  std::clog << "* filtration_simplex_range\n";
  for (auto sh : const_stree.filtration_simplex_range()) {
    std::clog << "   "
              << "[" << const_stree.filtration(sh) << "] ";
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* loop on simplex\n";
  for (std::size_t idx = 0; idx < const_stree.num_simplices(); idx++) {
    auto sh = const_stree.simplex(idx);
    std::clog << "   "
              << "[" << const_stree.filtration(sh) << "] ";
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* complex_simplex_range\n";
  for (auto sh : const_stree.complex_simplex_range()) {
    std::clog << "   "
              << "[" << const_stree.filtration(sh) << "] ";
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* complex_vertex_range\n";
  std::clog << "   ";
  for (auto vertex : const_stree.complex_vertex_range()) std::clog << "(" << vertex << ")";
  std::clog << std::endl;
  std::clog << "* skeleton_simplex_range in dim " << const_stree.dimension() << "\n";
  for (auto sh : const_stree.skeleton_simplex_range(const_stree.dimension())) {
    std::clog << "   "
              << "[" << const_stree.filtration(sh) << "] ";
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* star_simplex_range of {0, 1, 6}\n";
  for (auto sh : const_stree.star_simplex_range(const_stree.find({0, 1, 6}))) {
    std::clog << "   "
              << "[" << const_stree.filtration(sh) << "] ";
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* cofaces_simplex_range of {0, 1, 6} in codimension 1\n";
  for (auto sh : const_stree.cofaces_simplex_range(const_stree.find({0, 1, 6}), 1)) {
    std::clog << "   "
              << "[" << const_stree.filtration(sh) << "] ";
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* boundary_simplex_range of {0, 1, 6}\n";
  for (auto sh : const_stree.boundary_simplex_range(const_stree.find({0, 1, 6}))) {
    std::clog << "   "
              << "[" << const_stree.filtration(sh) << "] ";
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* boundary_opposite_vertex_simplex_range of {0, 1, 6}\n";
  for (auto face_opposite_vertex : const_stree.boundary_opposite_vertex_simplex_range(const_stree.find({0, 1, 6}))) {
    auto sh = face_opposite_vertex.first;
    std::clog << "   "
              << "[" << const_stree.filtration(sh) << "] ";
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* Hasse diagram\n";
  const_stree.print_hasse(std::cout);
  std::clog << "* vertex_with_same_filtration of {5, 3}"
            << const_stree.vertex_with_same_filtration(const_stree.find({5, 3})) << "\n";
  std::clog << "* edge_with_same_filtration of {0, 1, 2}\n";
  {
    auto sh = const_stree.edge_with_same_filtration(const_stree.find({0, 1, 2}));
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* minimal_simplex_with_same_filtration of {0, 1, 2}\n";
  {
    auto sh = const_stree.minimal_simplex_with_same_filtration(const_stree.find({0, 1, 2}));
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
  }
  std::clog << "* serialization size = " << const_stree.get_serialization_size() << "\n";
  std::clog << "* operator<<\n";
  std::clog << const_stree;

  std::clog << "* for_each_simplex\n";
  const_stree.for_each_simplex([&](auto sh, int dim) {
    std::clog << "dim = " << dim << " - ";
    for (auto vertex : const_stree.simplex_vertex_range(sh)) std::clog << "(" << vertex << ")";
    std::clog << std::endl;
    });

  BOOST_CHECK(const_stree.contiguous_vertices());
  BOOST_CHECK(const_stree.has_children(const_stree.find({0, 1})));
  [[maybe_unused]] auto sh = const_stree.root();                                    // -> OK
}

BOOST_AUTO_TEST_CASE(const_simplex_tree) {
  Simplex_tree<> st;

  st.insert_simplex_and_subfaces({2, 1, 0}, 3.0);
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, 4.0);
  st.insert_simplex_and_subfaces({3, 0}, 2.0);
  st.insert_simplex_and_subfaces({3, 4, 5}, 3.0);
  st.insert_simplex_and_subfaces({8}, 1.0);
  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o   o       */
  /*  2   0   3\X/4   8       */
  /*            o             */
  /*            5             */
  /*                          */

  // Required before browsing through filtration values
  st.initialize_filtration();

  test_simplex_tree_constness(st);

}

template<typename Simplex_tree>
void test_simplex_tree_data_constness(const Simplex_tree& const_stree) {
  std::clog << "* test_simplex_tree_data_constness\n";
  std::clog << "simplex_data({0, 1}) = " << const_stree.simplex_data(const_stree.find({0, 1})) << "\n";
  BOOST_CHECK(const_stree.simplex_data(const_stree.find({0, 1})) == std::string("{0, 1}"));
  std::clog << "simplex_data({0, 1, 2}) = " << const_stree.simplex_data(const_stree.find({0, 1, 2})) << "\n";
  BOOST_CHECK(const_stree.simplex_data(const_stree.find({2, 1, 0})) == std::string("{0, 1, 2}"));
}

struct Options_with_int_data : Simplex_tree_options_minimal {
  typedef std::string Simplex_data;
};

BOOST_AUTO_TEST_CASE(const_simplex_data) {
  Simplex_tree<Options_with_int_data> st;
  st.insert_simplex_and_subfaces({0, 1});
  st.insert_simplex_and_subfaces({2, 1});
  st.insert_simplex_and_subfaces({0, 2});
  st.simplex_data(st.find({0, 1})) = std::string("{0, 1}");
  st.expansion(3);
  st.simplex_data(st.find({0, 1, 2})) = std::string("{0, 1, 2}");
  test_simplex_tree_data_constness(st);
}
