#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_iostream_operator"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

struct MyOptions : Simplex_tree_options_full_featured {
  // Not doing persistence, so we don't need those
  static const bool store_key = false;
  static const bool store_filtration = false;
  // I have few vertices
  typedef short Vertex_handle;
};

typedef boost::mpl::list<Simplex_tree<>,
                         Simplex_tree<Simplex_tree_options_fast_persistence>
                        > list_of_tested_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(iostream_operator, Stree_type, list_of_tested_variants) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "SIMPLEX TREE IOSTREAM OPERATOR" << std::endl;

  Stree_type st;

  st.insert_simplex_and_subfaces({0, 1, 6, 7}, 4.0);
  st.insert_simplex_and_subfaces({3, 4, 5}, 3.0);
  st.insert_simplex_and_subfaces({3, 0}, 2.0);
  st.insert_simplex_and_subfaces({2, 1, 0}, 3.0);
  // FIXME
  st.set_filtration(4.0);

  st.initialize_filtration();
  // Display the Simplex_tree
  std::cout << "The ORIGINAL complex contains " << st.num_simplices() << " simplices" << std::endl;
  std::cout << "   - dimension " << st.dimension() << "   - filtration " << st.filtration() << std::endl;
  std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::cout << (int) vertex << " ";
    }
    std::cout << std::endl;
  }

  // st:
  //    1   6
  //    o---o
  //   /X\7/
  //  o---o---o---o
  //  2   0   3\X/4
  //            o
  //            5
  std::string iostream_file("simplex_tree_for_iostream_operator_unit_test.txt");
  std::ofstream simplex_tree_ostream(iostream_file.c_str());
  simplex_tree_ostream << st;
  simplex_tree_ostream.close();

  Stree_type read_st;
  std::ifstream simplex_tree_istream(iostream_file.c_str());
  simplex_tree_istream >> read_st;

  // Display the Simplex_tree
  std::cout << "The READ complex contains " << read_st.num_simplices() << " simplices" << std::endl;
  std::cout << "   - dimension " << read_st.dimension() << "   - filtration " << read_st.filtration() << std::endl;
  std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : read_st.filtration_simplex_range()) {
    std::cout << "   " << "[" << read_st.filtration(f_simplex) << "] ";
    for (auto vertex : read_st.simplex_vertex_range(f_simplex)) {
      std::cout << (int) vertex << " ";
    }
    std::cout << std::endl;
  }

  BOOST_CHECK(st == read_st);
}


BOOST_AUTO_TEST_CASE(mini_iostream_operator) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "MINI SIMPLEX TREE IOSTREAM OPERATOR" << std::endl;

  Simplex_tree<MyOptions> st;

  st.insert_simplex_and_subfaces({0, 1, 6, 7});
  st.insert_simplex_and_subfaces({3, 4, 5});
  st.insert_simplex_and_subfaces({3, 0});
  st.insert_simplex_and_subfaces({2, 1, 0});

  st.initialize_filtration();
  // Display the Simplex_tree
  std::cout << "The ORIGINAL complex contains " << st.num_simplices() << " simplices" << std::endl;
  std::cout << "   - dimension " << st.dimension() << "   - filtration " << st.filtration() << std::endl;
  std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::cout << (int) vertex << " ";
    }
    std::cout << std::endl;
  }

  // st:
  //    1   6
  //    o---o
  //   /X\7/
  //  o---o---o---o
  //  2   0   3\X/4
  //            o
  //            5
  std::string iostream_file("simplex_tree_for_iostream_operator_unit_test.txt");
  std::ofstream simplex_tree_ostream(iostream_file.c_str());
  simplex_tree_ostream << st;
  simplex_tree_ostream.close();

  Simplex_tree<MyOptions> read_st;
  std::ifstream simplex_tree_istream(iostream_file.c_str());
  simplex_tree_istream >> read_st;

  // Display the Simplex_tree
  std::cout << "The READ complex contains " << read_st.num_simplices() << " simplices" << std::endl;
  std::cout << "   - dimension " << read_st.dimension() << "   - filtration " << read_st.filtration() << std::endl;
  std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : read_st.filtration_simplex_range()) {
    std::cout << "   " << "[" << read_st.filtration(f_simplex) << "] ";
    for (auto vertex : read_st.simplex_vertex_range(f_simplex)) {
      std::cout << (int) vertex << " ";
    }
    std::cout << std::endl;
  }

  BOOST_CHECK(st == read_st);
}
