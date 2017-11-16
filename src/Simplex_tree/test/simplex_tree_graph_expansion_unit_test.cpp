#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility> // std::pair, std::make_pair
#include <cmath> // float comparison
#include <limits>
#include <functional> // greater

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

typedef boost::mpl::list<Simplex_tree<>, Simplex_tree<Simplex_tree_options_fast_persistence>> list_of_tested_variants;


bool AreAlmostTheSame(float a, float b) {
  return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion_with_blockers_3, typeST, list_of_tested_variants) {
  using Simplex_handle = typename typeST::Simplex_handle;
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  simplex_tree.expansion_with_blockers(3, [&](Simplex_handle sh){
      bool result = false;
      std::cout << "Blocker on [";
      // User can loop on the vertices from the given simplex_handle i.e.
      for (auto vertex : simplex_tree.simplex_vertex_range(sh)) {
        // We block the expansion, if the vertex '6' is in the given list of vertices
        if (vertex == 6)
          result = true;
        std::cout << vertex << ", ";
      }
      std::cout << "] ( " << simplex_tree.filtration(sh);
      // User can re-assign a new filtration value directly in the blocker (default is the maximal value of boudaries)
      simplex_tree.assign_filtration(sh, simplex_tree.filtration(sh) + 1.);

      std::cout << " + 1. ) = " << result << std::endl;

      return result;
    });

  std::cout << "********************************************************************\n";
  std::cout << "simplex_tree_expansion_with_blockers_3\n";
  std::cout << "********************************************************************\n";
  std::cout << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::cout << " - dimension " << simplex_tree.dimension() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::cout << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::cout << "(" << vertex << ")";
    std::cout << std::endl;
  }

  BOOST_CHECK(simplex_tree.num_simplices() == 23);
  BOOST_CHECK(simplex_tree.dimension() == 3);
  // {4, 5, 6} shall be blocked
  BOOST_CHECK(simplex_tree.find({4, 5, 6}) == simplex_tree.null_simplex());
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,2})), 4.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,3})), 5.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,2,3})), 6.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({1,2,3})), 6.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,2,3})), 7.));

}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion_with_blockers_2, typeST, list_of_tested_variants) {
  using Simplex_handle = typename typeST::Simplex_handle;
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  simplex_tree.expansion_with_blockers(2, [&](Simplex_handle sh){
      bool result = false;
      std::cout << "Blocker on [";
      // User can loop on the vertices from the given simplex_handle i.e.
      for (auto vertex : simplex_tree.simplex_vertex_range(sh)) {
        // We block the expansion, if the vertex '6' is in the given list of vertices
        if (vertex == 6)
          result = true;
        std::cout << vertex << ", ";
      }
      std::cout << "] ( " << simplex_tree.filtration(sh);
      // User can re-assign a new filtration value directly in the blocker (default is the maximal value of boudaries)
      simplex_tree.assign_filtration(sh, simplex_tree.filtration(sh) + 1.);

      std::cout << " + 1. ) = " << result << std::endl;

      return result;
    });

  std::cout << "********************************************************************\n";
  std::cout << "simplex_tree_expansion_with_blockers_2\n";
  std::cout << "********************************************************************\n";
  std::cout << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::cout << " - dimension " << simplex_tree.dimension() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::cout << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::cout << "(" << vertex << ")";
    std::cout << std::endl;
  }

  BOOST_CHECK(simplex_tree.num_simplices() == 22);
  BOOST_CHECK(simplex_tree.dimension() == 2);
  // {4, 5, 6} shall be blocked
  BOOST_CHECK(simplex_tree.find({4, 5, 6}) == simplex_tree.null_simplex());
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,2})), 4.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,3})), 5.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,2,3})), 6.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({1,2,3})), 6.));
  BOOST_CHECK(simplex_tree.find({0,1,2,3}) == simplex_tree.null_simplex());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion, typeST, list_of_tested_variants) {
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  simplex_tree.expansion(3);
  std::cout << "********************************************************************\n";
  std::cout << "simplex_tree_expansion_3\n";
  std::cout << "********************************************************************\n";
  std::cout << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::cout << " - dimension " << simplex_tree.dimension() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::cout << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::cout << "(" << vertex << ")";
    std::cout << std::endl;
  }

  BOOST_CHECK(simplex_tree.num_simplices() == 24);
  BOOST_CHECK(simplex_tree.dimension() == 3);

  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({4,5,6})), 10.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,2})), 3.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,3})), 4.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,2,3})), 5.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({1,2,3})), 5.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,2,3})), 5.));

}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_expansion_2, typeST, list_of_tested_variants) {
  // Construct the Simplex Tree with a 1-skeleton graph example
  typeST simplex_tree;

  simplex_tree.insert_simplex({0, 1}, 0.);
  simplex_tree.insert_simplex({0, 2}, 1.);
  simplex_tree.insert_simplex({0, 3}, 2.);
  simplex_tree.insert_simplex({1, 2}, 3.);
  simplex_tree.insert_simplex({1, 3}, 4.);
  simplex_tree.insert_simplex({2, 3}, 5.);
  simplex_tree.insert_simplex({2, 4}, 6.);
  simplex_tree.insert_simplex({3, 6}, 7.);
  simplex_tree.insert_simplex({4, 5}, 8.);
  simplex_tree.insert_simplex({4, 6}, 9.);
  simplex_tree.insert_simplex({5, 6}, 10.);
  simplex_tree.insert_simplex({6}, 10.);

  simplex_tree.expansion(2);

  std::cout << "********************************************************************\n";
  std::cout << "simplex_tree_expansion_2\n";
  std::cout << "********************************************************************\n";
  std::cout << "* The complex contains " << simplex_tree.num_simplices() << " simplices";
  std::cout << " - dimension " << simplex_tree.dimension() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::cout << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex))
      std::cout << "(" << vertex << ")";
    std::cout << std::endl;
  }

  BOOST_CHECK(simplex_tree.num_simplices() == 23);
  BOOST_CHECK(simplex_tree.dimension() == 2);

  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({4,5,6})), 10.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,2})), 3.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,1,3})), 4.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({0,2,3})), 5.));
  BOOST_CHECK(AreAlmostTheSame(simplex_tree.filtration(simplex_tree.find({1,2,3})), 5.));
  BOOST_CHECK(simplex_tree.find({0,1,2,3}) == simplex_tree.null_simplex());
}
