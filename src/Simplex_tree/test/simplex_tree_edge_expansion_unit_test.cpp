/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <vector>
#include <random>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_edge_expansion"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "gudhi/Simplex_tree.h"
#include <gudhi/Unitary_tests_utils.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

using namespace Gudhi;

struct Simplex_tree_options_fast_cofaces {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  static const bool link_nodes_by_label = true;
  static const bool stable_simplex_handles = false;
};

struct Simplex_tree_options_stable_simplex_handles_fast_cofaces : Simplex_tree_options_fast_cofaces {
  static const bool stable_simplex_handles = true;
};

using Point = std::vector<double>;

std::vector<Point> build_point_cloud(unsigned int numberOfPoints, int seed){
  std::vector<Point> finalPoints;
  std::set<Point> points;
  std::random_device dev;
  std::mt19937 rng(dev());
  if (seed > -1) rng.seed(seed);
  std::uniform_real_distribution<double> dist(0,10);

  for (unsigned int i = 0; i < numberOfPoints; ++i){
    auto res = points.insert({dist(rng), dist(rng)});
    while(!res.second){
      res = points.insert({dist(rng), dist(rng)});
    }
    finalPoints.push_back(*res.first);
  }

  return finalPoints;
}

template<class ST_type, typename St_options>
void test_insert_as_flag(ST_type& simplex_tree, int maxDim = -1) {
  Simplex_tree<St_options> st_to_test;
  std::vector<typename Simplex_tree<St_options>::Simplex_handle> added_simplices;

  for (auto& sh : simplex_tree.filtration_simplex_range()) {
    if (simplex_tree.dimension(sh) == 0) {
      auto u = *simplex_tree.simplex_vertex_range(sh).begin();
      st_to_test.insert_edge_as_flag(u, u, simplex_tree.filtration(sh), maxDim, added_simplices);
    } else if (simplex_tree.dimension(sh) == 1) {
      auto it = simplex_tree.simplex_vertex_range(sh).begin();
      auto u = *it;
      auto v = *(++it);
      st_to_test.insert_edge_as_flag(u, v, simplex_tree.filtration(sh), maxDim, added_simplices);
    }
  }

  BOOST_CHECK(added_simplices.size() == simplex_tree.num_simplices());
  BOOST_CHECK(st_to_test.dimension() == simplex_tree.dimension());
  BOOST_CHECK(st_to_test == simplex_tree);
}

template <class ST_type, typename St_options>
void test_unordered_insert_as_flag(ST_type& simplex_tree, int maxDim = -1) {
  Simplex_tree<St_options> st_to_test;
  std::vector<typename Simplex_tree<St_options>::Simplex_handle> added_simplices;

  // complex_simplex_range gives a lexicographical order where the word "12" is shorter than "1",
  // so vertices have to be inserted separately
  for (auto& sh : simplex_tree.skeleton_simplex_range(0)) {
    auto u = *simplex_tree.simplex_vertex_range(sh).begin();
    st_to_test.insert_edge_as_flag(u, u, simplex_tree.filtration(sh), maxDim, added_simplices);
  }

  for (auto& sh : simplex_tree.complex_simplex_range()) {
    if (simplex_tree.dimension(sh) == 1) {
      auto it = simplex_tree.simplex_vertex_range(sh).begin();
      auto u = *it;
      auto v = *(++it);
      st_to_test.insert_edge_as_flag(u, v, simplex_tree.filtration(sh), maxDim, added_simplices);
    }
  }
  // as the insertions are not in the order of the filtration,
  // the filtration values for higher dimensional simplices has to be restored
  st_to_test.make_filtration_non_decreasing();

  BOOST_CHECK(added_simplices.size() == simplex_tree.num_simplices());
  BOOST_CHECK(st_to_test.dimension() == simplex_tree.dimension());
  BOOST_CHECK(st_to_test == simplex_tree);
}

BOOST_AUTO_TEST_CASE(simplex_tree_random_rips_expansion) {
  using Filtration_value = typename Simplex_tree<>::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

  for (unsigned int i = 10; i < 100; i += 5){
    Simplex_tree<> simplex_tree;
    Rips_complex rips_complex(build_point_cloud(i, -1), 3, Gudhi::Euclidean_distance());
    rips_complex.create_complex(simplex_tree, 100);

    std::clog << i << "  ********************************************************************\n";
    std::clog << "simplex_tree_random_rips_expansion\n";
    std::clog << "Test tree 1: insertion by filtration values\n";
    std::clog << "************************************************************************\n";
    test_insert_as_flag<Simplex_tree<>,Simplex_tree_options_stable_simplex_handles_fast_cofaces>(simplex_tree);
    test_insert_as_flag<Simplex_tree<>,Simplex_tree_options_fast_cofaces>(simplex_tree);

    std::clog << i << "  ********************************************************************\n";
    std::clog << "simplex_tree_random_rips_expansion\n";
    std::clog << "Test tree 2: insertion by lexicographical order\n"; //equivalent to random order
    std::clog << "************************************************************************\n";
    test_unordered_insert_as_flag<Simplex_tree<>,Simplex_tree_options_stable_simplex_handles_fast_cofaces>(simplex_tree);
    test_unordered_insert_as_flag<Simplex_tree<>,Simplex_tree_options_fast_cofaces>(simplex_tree);
  }
}

BOOST_AUTO_TEST_CASE(simplex_tree_random_rips_expansion_with_max_dim) {
  using Filtration_value = typename Simplex_tree<>::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

  for (unsigned int i = 10; i < 100; i += 5){
    int maxDim = i / 5;
    Simplex_tree<> simplex_tree;
    Rips_complex rips_complex(build_point_cloud(i, -1), 3, Gudhi::Euclidean_distance());
    rips_complex.create_complex(simplex_tree, maxDim);

    std::clog << i << "  ********************************************************************\n";
    std::clog << "simplex_tree_random_rips_expansion\n";
    std::clog << "Test tree 3: insertion by filtration values and max dimension = " << maxDim << "\n";
    std::clog << "************************************************************************\n";
    test_insert_as_flag<Simplex_tree<>,Simplex_tree_options_stable_simplex_handles_fast_cofaces>(simplex_tree, maxDim);
    test_insert_as_flag<Simplex_tree<>,Simplex_tree_options_fast_cofaces>(simplex_tree, maxDim);

    std::clog << i << "  ********************************************************************\n";
    std::clog << "simplex_tree_random_rips_expansion\n";
    std::clog << "Test tree 4: insertion by lexicographical order and max dimension = " << maxDim << "\n";   //equivalent to random order
    std::clog << "************************************************************************\n";
    test_unordered_insert_as_flag<Simplex_tree<>,Simplex_tree_options_stable_simplex_handles_fast_cofaces>(simplex_tree, maxDim);
    test_unordered_insert_as_flag<Simplex_tree<>,Simplex_tree_options_fast_cofaces>(simplex_tree, maxDim);
  }
}

BOOST_AUTO_TEST_CASE(flag_expansion) {
  using typeST = Simplex_tree<Simplex_tree_options_stable_simplex_handles_fast_cofaces>;
  {
    std::cout << "************************************************************************" << std::endl;
    std::cout << "Test flag expansion 1" << std::endl;
    std::cout << "************************************************************************" << std::endl;

    typeST st;
    std::vector<typename typeST::Simplex_handle> added_simplices;

    st.insert_edge_as_flag(0,0,0,6,added_simplices);
    st.insert_edge_as_flag(1,1,0,6,added_simplices);
    st.insert_edge_as_flag(2,2,0,6,added_simplices);
    st.insert_edge_as_flag(3,3,0,6,added_simplices);
    st.insert_edge_as_flag(4,4,0,6,added_simplices);
    st.insert_edge_as_flag(5,5,0,6,added_simplices);
    st.insert_edge_as_flag(6,6,0,6,added_simplices);
    st.insert_edge_as_flag(0,1,0,6,added_simplices);
    st.insert_edge_as_flag(0,2,0,6,added_simplices);
    st.insert_edge_as_flag(0,3,0,6,added_simplices);
    st.insert_edge_as_flag(0,4,0,6,added_simplices);
    st.insert_edge_as_flag(0,5,0,6,added_simplices);
    st.insert_edge_as_flag(0,6,0,6,added_simplices);
    st.insert_edge_as_flag(1,2,0,6,added_simplices);
    st.insert_edge_as_flag(1,3,0,6,added_simplices);
    st.insert_edge_as_flag(1,4,0,6,added_simplices);
    st.insert_edge_as_flag(1,5,0,6,added_simplices);
    st.insert_edge_as_flag(1,6,0,6,added_simplices);
    st.insert_edge_as_flag(2,3,0,6,added_simplices);
    st.insert_edge_as_flag(2,4,0,6,added_simplices);
    st.insert_edge_as_flag(2,5,0,6,added_simplices);
    st.insert_edge_as_flag(2,6,0,6,added_simplices);
    st.insert_edge_as_flag(3,4,0,6,added_simplices);
    st.insert_edge_as_flag(3,5,0,6,added_simplices);
    st.insert_edge_as_flag(3,6,0,6,added_simplices);
    st.insert_edge_as_flag(4,5,0,6,added_simplices);
    st.insert_edge_as_flag(4,6,0,6,added_simplices);
    st.insert_edge_as_flag(5,6,0,6,added_simplices);

    st.assign_filtration(st.find({0,2,4}), 10);
    st.assign_filtration(st.find({1,5}), 20);
    st.assign_filtration(st.find({1,2,4}), 30);
    st.assign_filtration(st.find({3}), 5);
    st.make_filtration_non_decreasing();

    BOOST_CHECK_EQUAL(added_simplices.size(), 127);
    BOOST_CHECK(st.filtration(st.find({1,2}))==0);
    BOOST_CHECK(st.filtration(st.find({0,1,2,3,4}))==30);
    BOOST_CHECK(st.minimal_simplex_with_same_filtration(st.find({0,1,2,3,4,5}))==st.find({1,2,4}));
    BOOST_CHECK(st.minimal_simplex_with_same_filtration(st.find({0,2,3}))==st.find({3}));
    auto s=st.minimal_simplex_with_same_filtration(st.find({0,2,6}));
    BOOST_CHECK(s==st.find({0})||s==st.find({2})||s==st.find({6}));
    BOOST_CHECK(st.vertex_with_same_filtration(st.find({2}))==2);
    BOOST_CHECK(st.vertex_with_same_filtration(st.find({1,5}))==st.null_vertex());
    BOOST_CHECK(st.vertex_with_same_filtration(st.find({5,6}))>=5);
  }
  {
    std::cout << "************************************************************************" << std::endl;
    std::cout << "Test flag expansion 2" << std::endl;
    std::cout << "************************************************************************" << std::endl;

    typeST st;
    std::vector<typename typeST::Simplex_handle> added_simplices;

    st.insert_edge_as_flag(0,0,0,6,added_simplices);
    st.insert_edge_as_flag(1,1,0,6,added_simplices);
    st.insert_edge_as_flag(2,2,0,6,added_simplices);
    st.insert_edge_as_flag(4,4,0,6,added_simplices);
    st.insert_edge_as_flag(5,5,0,6,added_simplices);
    st.insert_edge_as_flag(6,6,0,6,added_simplices);

    st.insert_edge_as_flag(0,1,0,6,added_simplices);
    st.insert_edge_as_flag(0,5,0,6,added_simplices);
    st.insert_edge_as_flag(0,6,0,6,added_simplices);
    st.insert_edge_as_flag(1,6,0,6,added_simplices);
    st.insert_edge_as_flag(2,5,0,6,added_simplices);
    st.insert_edge_as_flag(2,6,0,6,added_simplices);
    st.insert_edge_as_flag(4,5,0,6,added_simplices);
    st.insert_edge_as_flag(4,6,0,6,added_simplices);
    st.insert_edge_as_flag(5,6,0,6,added_simplices);

    st.insert_edge_as_flag(3,3,5,6,added_simplices);
    st.insert_edge_as_flag(0,3,5,6,added_simplices);
    st.insert_edge_as_flag(1,3,5,6,added_simplices);
    st.insert_edge_as_flag(2,3,5,6,added_simplices);
    st.insert_edge_as_flag(3,4,5,6,added_simplices);
    st.insert_edge_as_flag(3,5,5,6,added_simplices);
    st.insert_edge_as_flag(3,6,5,6,added_simplices);

    st.insert_edge_as_flag(0,2,10,6,added_simplices);
    st.insert_edge_as_flag(0,4,10,6,added_simplices);
    st.insert_edge_as_flag(1,5,20,6,added_simplices);
    st.insert_edge_as_flag(1,2,30,6,added_simplices);
    st.insert_edge_as_flag(1,4,30,6,added_simplices);
    st.insert_edge_as_flag(2,4,30,6,added_simplices);

    BOOST_CHECK_EQUAL(added_simplices.size(), 127);
    BOOST_CHECK(st.filtration(st.find({1,2}))==30);
    BOOST_CHECK(st.filtration(st.find({0,1,2,3,4}))==30);
    BOOST_CHECK(st.minimal_simplex_with_same_filtration(st.find({0,1,2,3,4,5}))==st.find({1,2}));
    BOOST_CHECK(st.minimal_simplex_with_same_filtration(st.find({0,2,3}))==st.find({0,2}));
    auto s=st.minimal_simplex_with_same_filtration(st.find({0,2,6}));
    BOOST_CHECK(s==st.find({0,2}));
    BOOST_CHECK(st.vertex_with_same_filtration(st.find({2}))==2);
    BOOST_CHECK(st.vertex_with_same_filtration(st.find({1,5}))==st.null_vertex());
    BOOST_CHECK(st.vertex_with_same_filtration(st.find({5,6}))>=5);
  }
  {
    std::cout << "************************************************************************" << std::endl;
    std::cout << "Test flag expansion 3" << std::endl;
    std::cout << "************************************************************************" << std::endl;

    typeST st;
    std::vector<typename typeST::Simplex_handle> added_simplices;

    st.insert_edge_as_flag(0,0,0,6,added_simplices);

    st.insert_edge_as_flag(1,1,0,6,added_simplices);
    st.insert_edge_as_flag(0,1,0,6,added_simplices);

    st.insert_edge_as_flag(2,2,0,6,added_simplices);
    st.insert_edge_as_flag(0,2,10,6,added_simplices);
    st.insert_edge_as_flag(1,2,30,6,added_simplices);

    st.insert_edge_as_flag(3,3,5,6,added_simplices);
    st.insert_edge_as_flag(0,3,5,6,added_simplices);
    st.insert_edge_as_flag(1,3,5,6,added_simplices);
    st.insert_edge_as_flag(2,3,5,6,added_simplices);

    st.insert_edge_as_flag(4,4,0,6,added_simplices);
    st.insert_edge_as_flag(0,4,10,6,added_simplices);
    st.insert_edge_as_flag(1,4,30,6,added_simplices);
    st.insert_edge_as_flag(2,4,30,6,added_simplices);
    st.insert_edge_as_flag(3,4,5,6,added_simplices);

    st.insert_edge_as_flag(5,5,0,6,added_simplices);
    st.insert_edge_as_flag(0,5,0,6,added_simplices);
    st.insert_edge_as_flag(1,5,20,6,added_simplices);
    st.insert_edge_as_flag(2,5,0,6,added_simplices);
    st.insert_edge_as_flag(3,5,5,6,added_simplices);
    st.insert_edge_as_flag(4,5,0,6,added_simplices);

    st.insert_edge_as_flag(6,6,0,6,added_simplices);
    st.insert_edge_as_flag(0,6,0,6,added_simplices);
    st.insert_edge_as_flag(1,6,0,6,added_simplices);
    st.insert_edge_as_flag(2,6,0,6,added_simplices);
    st.insert_edge_as_flag(3,6,5,6,added_simplices);
    st.insert_edge_as_flag(4,6,0,6,added_simplices);
    st.insert_edge_as_flag(5,6,0,6,added_simplices);

    st.make_filtration_non_decreasing();

    BOOST_CHECK_EQUAL(added_simplices.size(), 127);
    BOOST_CHECK(st.filtration(st.find({1,2}))==30);
    BOOST_CHECK(st.filtration(st.find({0,1,2,3,4}))==30);
    BOOST_CHECK(st.minimal_simplex_with_same_filtration(st.find({0,1,2,3,4,5}))==st.find({1,2}));
    BOOST_CHECK(st.minimal_simplex_with_same_filtration(st.find({0,2,3}))==st.find({0,2}));
    auto s=st.minimal_simplex_with_same_filtration(st.find({0,2,6}));
    BOOST_CHECK(s==st.find({0,2}));
    BOOST_CHECK(st.vertex_with_same_filtration(st.find({2}))==2);
    BOOST_CHECK(st.vertex_with_same_filtration(st.find({1,5}))==st.null_vertex());
    BOOST_CHECK(st.vertex_with_same_filtration(st.find({5,6}))>=5);
  }
  {
    std::cout << "************************************************************************" << std::endl;
    std::cout << "Test flag expansion 4" << std::endl;
    std::cout << "************************************************************************" << std::endl;

    typeST st;
    std::vector<typename typeST::Simplex_handle> added_simplices;

    st.insert_edge_as_flag(2,2,2,50,added_simplices);
    st.insert_edge_as_flag(5,5,2,50,added_simplices);
    st.insert_edge_as_flag(2,5,2,50,added_simplices);

    st.insert_edge_as_flag(0,0,3,50,added_simplices);
    st.insert_edge_as_flag(0,5,3,50,added_simplices);

    st.insert_edge_as_flag(1,1,4,50,added_simplices);
    st.insert_edge_as_flag(1,5,4,50,added_simplices);

    st.insert_edge_as_flag(1,2,5,50,added_simplices);

    st.insert_edge_as_flag(3,3,6,50,added_simplices);
    st.insert_edge_as_flag(4,4,6,50,added_simplices);
    st.insert_edge_as_flag(3,4,6,50,added_simplices);

    st.insert_edge_as_flag(0,1,8,50,added_simplices);

    st.insert_edge_as_flag(1,3,9,50,added_simplices);

    st.insert_edge_as_flag(0,2,10,50,added_simplices);

    BOOST_CHECK_EQUAL(added_simplices.size(), 19);
    BOOST_CHECK(st.edge_with_same_filtration(st.find({0,1,2,5}))==st.find({0,2}));
    BOOST_CHECK(st.edge_with_same_filtration(st.find({1,5}))==st.find({1,5}));
  }
}
