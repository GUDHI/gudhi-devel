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

typedef boost::mpl::list<Simplex_tree<>,
                         Simplex_tree<Simplex_tree_options_fast_persistence>,
                         Simplex_tree<Simplex_tree_options_fast_cofaces>> list_of_tested_variants;

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

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_random_rips_expansion1, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_random_rips_expansion\n";
  std::clog << "Test tree 1: insertion by filtration values\n";
  std::clog << "********************************************************************\n";

  using Filtration_value = typename typeST::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

  for (unsigned int i = 10; i < 100; i += 5){
    typeST simplex_tree;
    Rips_complex rips_complex(build_point_cloud(i, -1), 3, Gudhi::Euclidean_distance());
    rips_complex.create_complex(simplex_tree, 100);

    Simplex_tree<Simplex_tree_options_fast_cofaces> st_to_test;
    std::vector<Simplex_tree<Simplex_tree_options_fast_cofaces>::Simplex_handle> added_simplices;

    for (auto& sh : simplex_tree.filtration_simplex_range()){
      if (simplex_tree.dimension(sh) == 0){
        auto u = *simplex_tree.simplex_vertex_range(sh).begin();
        st_to_test.insert_edge_as_flag(u, u, simplex_tree.filtration(sh), -1, added_simplices);
      } else if (simplex_tree.dimension(sh) == 1){
        auto it = simplex_tree.simplex_vertex_range(sh).begin();
        auto u = *it;
        auto v = *(++it);
        st_to_test.insert_edge_as_flag(u, v, simplex_tree.filtration(sh), -1, added_simplices);
      }
    }

    BOOST_CHECK(added_simplices.size() == simplex_tree.num_simplices());
    BOOST_CHECK(st_to_test == simplex_tree);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_random_rips_expansion2, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_random_rips_expansion\n";
  std::clog << "Test tree 2: insertion by lexicographical order\n";
  std::clog << "********************************************************************\n";

  using Filtration_value = typename typeST::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

  for (unsigned int i = 10; i < 100; i += 5){
    typeST simplex_tree;
    Rips_complex rips_complex(build_point_cloud(i, -1), 3, Gudhi::Euclidean_distance());
    rips_complex.create_complex(simplex_tree, 100);

    Simplex_tree<Simplex_tree_options_fast_cofaces> st_to_test;
    std::vector<Simplex_tree<Simplex_tree_options_fast_cofaces>::Simplex_handle> added_simplices;

    // complex_simplex_range gives a lexicographical order where the word "12" is shorter than "1",
    // so vertices have to be inserted separately
    for (auto& sh : simplex_tree.skeleton_simplex_range(0)){
      auto u = *simplex_tree.simplex_vertex_range(sh).begin();
      st_to_test.insert_edge_as_flag(u, u, simplex_tree.filtration(sh), -1, added_simplices);
    }

    for (auto& sh : simplex_tree.complex_simplex_range()){
      if (simplex_tree.dimension(sh) == 1){
        auto it = simplex_tree.simplex_vertex_range(sh).begin();
        auto u = *it;
        auto v = *(++it);
        st_to_test.insert_edge_as_flag(u, v, simplex_tree.filtration(sh), -1, added_simplices);
      }
    }
    // as the insertions are not in the order of the filtration,
    // the filtration values for higher dimensional simplices has to be restored
    st_to_test.make_filtration_non_decreasing();

    BOOST_CHECK(added_simplices.size() == simplex_tree.num_simplices());
    BOOST_CHECK(st_to_test == simplex_tree);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_random_rips_expansion1_with_max_dim, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_random_rips_expansion\n";
  std::clog << "Test tree 3: insertion by filtration values and max dimension\n";
  std::clog << "********************************************************************\n";

  using Filtration_value = typename typeST::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

  for (unsigned int i = 10; i < 100; i += 5){
    int maxDim = i / 5;
    typeST simplex_tree;
    Rips_complex rips_complex(build_point_cloud(i, -1), 3, Gudhi::Euclidean_distance());
    rips_complex.create_complex(simplex_tree, maxDim);

    Simplex_tree<Simplex_tree_options_fast_cofaces> st_to_test;
    std::vector<Simplex_tree<Simplex_tree_options_fast_cofaces>::Simplex_handle> added_simplices;

    for (auto& sh : simplex_tree.filtration_simplex_range()){
      if (simplex_tree.dimension(sh) == 0){
        auto u = *simplex_tree.simplex_vertex_range(sh).begin();
        st_to_test.insert_edge_as_flag(u, u, simplex_tree.filtration(sh), maxDim, added_simplices);
      } else if (simplex_tree.dimension(sh) == 1){
        auto it = simplex_tree.simplex_vertex_range(sh).begin();
        auto u = *it;
        auto v = *(++it);
        st_to_test.insert_edge_as_flag(u, v, simplex_tree.filtration(sh), maxDim, added_simplices);
      }
    }

    BOOST_CHECK(added_simplices.size() == simplex_tree.num_simplices());
    BOOST_CHECK(st_to_test == simplex_tree);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_random_rips_expansion2_with_max_dim, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************\n";
  std::clog << "simplex_tree_random_rips_expansion\n";
  std::clog << "Test tree 4: insertion by lexicographical order and max dimension\n";
  std::clog << "********************************************************************\n";

  using Filtration_value = typename typeST::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

  for (unsigned int i = 10; i < 100; i += 5){
    int maxDim = i / 5;
    typeST simplex_tree;
    Rips_complex rips_complex(build_point_cloud(i, -1), 3, Gudhi::Euclidean_distance());
    rips_complex.create_complex(simplex_tree, maxDim);

    Simplex_tree<Simplex_tree_options_fast_cofaces> st_to_test;
    std::vector<Simplex_tree<Simplex_tree_options_fast_cofaces>::Simplex_handle> added_simplices;

    // complex_simplex_range gives a lexicographical order where the word "12" is shorter than "1",
    // so vertices have to be inserted separately
    for (auto& sh : simplex_tree.skeleton_simplex_range(0)){
      auto u = *simplex_tree.simplex_vertex_range(sh).begin();
      st_to_test.insert_edge_as_flag(u, u, simplex_tree.filtration(sh), maxDim, added_simplices);
    }

    for (auto& sh : simplex_tree.complex_simplex_range()){
      if (simplex_tree.dimension(sh) == 1){
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
    BOOST_CHECK(st_to_test == simplex_tree);
  }
}





