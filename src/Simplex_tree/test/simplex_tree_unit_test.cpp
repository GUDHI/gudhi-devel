/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>  // std::pair, std::make_pair
#include <cmath>  // float comparison
#include <limits>
#include <functional>  // greater
#include <tuple>  // std::tie

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

typedef boost::mpl::list<Simplex_tree<>, Simplex_tree<Simplex_tree_options_fast_persistence>> list_of_tested_variants;


template<class typeST>
void test_empty_simplex_tree(typeST& tst) {
  typedef typename typeST::Vertex_handle Vertex_handle;
  const Vertex_handle DEFAULT_VERTEX_VALUE = Vertex_handle(- 1);
  BOOST_CHECK(tst.null_vertex() == DEFAULT_VERTEX_VALUE);
  BOOST_CHECK(tst.num_vertices() == (size_t) 0);
  BOOST_CHECK(tst.num_simplices() == (size_t) 0);
  typename typeST::Siblings* STRoot = tst.root();
  BOOST_CHECK(STRoot != nullptr);
  BOOST_CHECK(STRoot->oncles() == nullptr);
  BOOST_CHECK(STRoot->parent() == DEFAULT_VERTEX_VALUE);
  BOOST_CHECK(tst.dimension() == -1);
}

template<class typeST>
void test_iterators_on_empty_simplex_tree(typeST& tst) {
  std::cout << "Iterator on vertices: " << std::endl;
  for (auto vertex : tst.complex_vertex_range()) {
    std::cout << "vertice:" << vertex << std::endl;
    BOOST_CHECK(false); // shall be empty
  }
  std::cout << "Iterator on simplices: " << std::endl;
  for (auto simplex : tst.complex_simplex_range()) {
    BOOST_CHECK(simplex != simplex); // shall be empty - to remove warning of non-used simplex
  }

  std::cout
      << "Iterator on Simplices in the filtration, with [filtration value]:"
      << std::endl;
  for (auto f_simplex : tst.filtration_simplex_range()) {
    BOOST_CHECK(false); // shall be empty
    std::cout << "test_iterators_on_empty_simplex_tree - filtration="
        << tst.filtration(f_simplex) << std::endl;
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_when_empty, typeST, list_of_tested_variants) {
  typedef std::pair<typename typeST::Simplex_handle, bool> typePairSimplexBool;
  typedef std::vector<typename typeST::Vertex_handle> typeVectorVertex;

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF DEFAULT CONSTRUCTOR" << std::endl;
  typeST st;

  test_empty_simplex_tree(st);

  test_iterators_on_empty_simplex_tree(st);
  // TEST OF EMPTY INSERTION
  std::cout << "TEST OF EMPTY INSERTION" << std::endl;
  typeVectorVertex simplexVectorEmpty;
  BOOST_CHECK(simplexVectorEmpty.empty() == true);
  typePairSimplexBool returnEmptyValue = st.insert_simplex(simplexVectorEmpty, 0.0);
  BOOST_CHECK(returnEmptyValue.first == typename typeST::Simplex_handle(nullptr));
  BOOST_CHECK(returnEmptyValue.second == true);

  test_empty_simplex_tree(st);

  test_iterators_on_empty_simplex_tree(st);
}

bool AreAlmostTheSame(float a, float b) {
  return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_from_file, typeST, list_of_tested_variants) {
  // TEST OF INSERTION
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF SIMPLEX TREE FROM A FILE" << std::endl;
  typeST st;

  std::string inputFile("simplex_tree_for_unit_test.txt");
  std::ifstream simplex_tree_stream(inputFile.c_str());
  simplex_tree_stream >> st;

  // Display the Simplex_tree
  std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
  std::cout << "   - dimension " << st.dimension() << std::endl;

  // Check
  BOOST_CHECK(st.num_simplices() == 143353);
  BOOST_CHECK(st.dimension() == 3);

  int previous_size = 0;
  for (auto f_simplex : st.filtration_simplex_range()) {
    // Size of simplex
    int size = 0;
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      // Remove warning
      (void) vertex;
      size++;
    }
    BOOST_CHECK(AreAlmostTheSame(st.filtration(f_simplex), (0.1 * size))); // Specific test: filtration = 0.1 * simplex_size
    BOOST_CHECK(previous_size <= size); // Check list is sorted (because of sorted filtrations in simplex_tree.txt)
    previous_size = size;
  }
  simplex_tree_stream.close();
}

template<class typeST, class typeSimplex>
void test_simplex_tree_contains(typeST& simplexTree, typeSimplex& simplex, int pos) {
  auto f_simplex = simplexTree.filtration_simplex_range().begin() + pos;

  std::cout << "test_simplex_tree_contains - filtration=" << simplexTree.filtration(*f_simplex) << "||" << simplex.second << std::endl;
  BOOST_CHECK(AreAlmostTheSame(simplexTree.filtration(*f_simplex), simplex.second));

  int simplexIndex = simplex.first.size() - 1;
  std::sort(simplex.first.begin(), simplex.first.end()); // if the simplex wasn't sorted, the next test could fail
  for (auto vertex : simplexTree.simplex_vertex_range(*f_simplex)) {
    std::cout << "test_simplex_tree_contains - vertex=" << vertex << "||" << simplex.first.at(simplexIndex) << std::endl;
    BOOST_CHECK(vertex == simplex.first.at(simplexIndex));
    BOOST_CHECK(simplexIndex >= 0);
    simplexIndex--;
  }
}

template<class typeST, class typePairSimplexBool>
void test_simplex_tree_insert_returns_true(const typePairSimplexBool& returnValue) {
  BOOST_CHECK(returnValue.second == true);
  // Simplex_handle = boost::container::flat_map< typeST::Vertex_handle, Node >::iterator
  typename typeST::Simplex_handle shReturned = returnValue.first;
  BOOST_CHECK(shReturned != typename typeST::Simplex_handle(nullptr));
}

// Global variables
int dim_max = -1;

template<class typeST, class Filtration_value>
void set_and_test_simplex_tree_dim_fil(typeST& simplexTree, int vectorSize, const Filtration_value& fil) {
  if (vectorSize > dim_max + 1) {
    dim_max = vectorSize - 1;
    simplexTree.set_dimension(dim_max);
    std::cout << "   set_and_test_simplex_tree_dim_fil - dim_max=" << dim_max
        << std::endl;
  }

  BOOST_CHECK(simplexTree.dimension() == dim_max);

  // Another way to count simplices:
  size_t num_simp = 0;
  for (auto f_simplex : simplexTree.complex_simplex_range()) {
    // Remove warning
    (void) f_simplex;
    num_simp++;
  }

  BOOST_CHECK(simplexTree.num_simplices() == num_simp);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_insertion, typeST, list_of_tested_variants) {
  typedef typename typeST::Filtration_value Filtration_value;
  typedef std::pair<typename typeST::Simplex_handle, bool> typePairSimplexBool;
  typedef std::vector<typename typeST::Vertex_handle> typeVectorVertex;
  typedef std::pair<typeVectorVertex, Filtration_value> typeSimplex;
  const Filtration_value FIRST_FILTRATION_VALUE = 0.1;
  const Filtration_value SECOND_FILTRATION_VALUE = 0.2;
  const Filtration_value THIRD_FILTRATION_VALUE = 0.3;
  const Filtration_value FOURTH_FILTRATION_VALUE = 0.4;
  // reset since we run the test several times
  dim_max = -1;

  // TEST OF INSERTION
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF INSERTION" << std::endl;
  typeST st;

  // ++ FIRST
  std::cout << "   - INSERT 0" << std::endl;
  typeVectorVertex firstSimplexVector{0};
  BOOST_CHECK(firstSimplexVector.size() == 1);
  typeSimplex firstSimplex = std::make_pair(firstSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
  typePairSimplexBool returnValue = st.insert_simplex(firstSimplex.first, firstSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, firstSimplexVector.size(), firstSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 1);

  // ++ SECOND
  std::cout << "   - INSERT 1" << std::endl;
  typeVectorVertex secondSimplexVector{1};
  BOOST_CHECK(secondSimplexVector.size() == 1);
  typeSimplex secondSimplex = std::make_pair(secondSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
  returnValue = st.insert_simplex(secondSimplex.first, secondSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, secondSimplexVector.size(), secondSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 2);

  // ++ THIRD
  std::cout << "   - INSERT (0,1)" << std::endl;
  typeVectorVertex thirdSimplexVector{0, 1};
  BOOST_CHECK(thirdSimplexVector.size() == 2);
  typeSimplex thirdSimplex = std::make_pair(thirdSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
  returnValue = st.insert_simplex(thirdSimplex.first, thirdSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, thirdSimplexVector.size(), thirdSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 2); // Not incremented !!

  // ++ FOURTH
  std::cout << "   - INSERT 2" << std::endl;
  typeVectorVertex fourthSimplexVector{2};
  BOOST_CHECK(fourthSimplexVector.size() == 1);
  typeSimplex fourthSimplex = std::make_pair(fourthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
  returnValue = st.insert_simplex(fourthSimplex.first, fourthSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, fourthSimplexVector.size(), fourthSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 3);

  // ++ FIFTH
  std::cout << "   - INSERT (2,0)" << std::endl;
  typeVectorVertex fifthSimplexVector{2, 0};
  BOOST_CHECK(fifthSimplexVector.size() == 2);
  typeSimplex fifthSimplex = std::make_pair(fifthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
  returnValue = st.insert_simplex(fifthSimplex.first, fifthSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, fifthSimplexVector.size(), fifthSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 3); // Not incremented !!

  // ++ SIXTH
  std::cout << "   - INSERT (2,1)" << std::endl;
  typeVectorVertex sixthSimplexVector{2, 1};
  BOOST_CHECK(sixthSimplexVector.size() == 2);
  typeSimplex sixthSimplex = std::make_pair(sixthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
  returnValue = st.insert_simplex(sixthSimplex.first, sixthSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, sixthSimplexVector.size(), sixthSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 3); // Not incremented !!

  // ++ SEVENTH
  std::cout << "   - INSERT (2,1,0)" << std::endl;
  typeVectorVertex seventhSimplexVector{2, 1, 0};
  BOOST_CHECK(seventhSimplexVector.size() == 3);
  typeSimplex seventhSimplex = std::make_pair(seventhSimplexVector, Filtration_value(THIRD_FILTRATION_VALUE));
  returnValue = st.insert_simplex(seventhSimplex.first, seventhSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, seventhSimplexVector.size(), seventhSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 3); // Not incremented !!

  // ++ EIGHTH
  std::cout << "   - INSERT 3" << std::endl;
  typeVectorVertex eighthSimplexVector{3};
  BOOST_CHECK(eighthSimplexVector.size() == 1);
  typeSimplex eighthSimplex = std::make_pair(eighthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
  returnValue = st.insert_simplex(eighthSimplex.first, eighthSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, eighthSimplexVector.size(), eighthSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 4);

  // ++ NINETH
  std::cout << "   - INSERT (3,0)" << std::endl;
  typeVectorVertex ninethSimplexVector{3, 0};
  BOOST_CHECK(ninethSimplexVector.size() == 2);
  typeSimplex ninethSimplex = std::make_pair(ninethSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
  returnValue = st.insert_simplex(ninethSimplex.first, ninethSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, ninethSimplexVector.size(), ninethSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 4); // Not incremented !!

  // ++ TENTH
  std::cout << "   - INSERT 0 (already inserted)" << std::endl;
  typeVectorVertex tenthSimplexVector{0};
  BOOST_CHECK(tenthSimplexVector.size() == 1);
  // With a different filtration value
  typeSimplex tenthSimplex = std::make_pair(tenthSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE));
  returnValue = st.insert_simplex(tenthSimplex.first, tenthSimplex.second);

  BOOST_CHECK(returnValue.second == false);
  // Simplex_handle = boost::container::flat_map< typeST::Vertex_handle, Node >::iterator
  typename typeST::Simplex_handle shReturned = returnValue.first;
  BOOST_CHECK(shReturned == typename typeST::Simplex_handle(nullptr));
  std::cout << "st.num_vertices()=" << st.num_vertices() << std::endl;
  BOOST_CHECK(st.num_vertices() == (size_t) 4); // Not incremented !!
  BOOST_CHECK(st.dimension() == dim_max);

  // ++ ELEVENTH
  std::cout << "   - INSERT (2,1,0) (already inserted)" << std::endl;
  typeVectorVertex eleventhSimplexVector{2, 1, 0};
  BOOST_CHECK(eleventhSimplexVector.size() == 3);
  typeSimplex eleventhSimplex = std::make_pair(eleventhSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE));
  returnValue = st.insert_simplex(eleventhSimplex.first, eleventhSimplex.second);

  BOOST_CHECK(returnValue.second == false);
  // Simplex_handle = boost::container::flat_map< typeST::Vertex_handle, Node >::iterator
  shReturned = returnValue.first;
  BOOST_CHECK(shReturned == typename typeST::Simplex_handle(nullptr));
  BOOST_CHECK(st.num_vertices() == (size_t) 4); // Not incremented !!
  BOOST_CHECK(st.dimension() == dim_max);

  /* Inserted simplex:        */
  /*    1                     */
  /*    o                     */
  /*   /X\                    */
  /*  o---o---o               */
  /*  2   0   3               */

  //   [0.1] 0
  //   [0.1] 1
  //   [0.1] 2
  //   [0.1] 3
  //   [0.2] 1 0
  //   [0.2] 2 0
  //   [0.2] 2 1
  //   [0.2] 3 0
  //   [0.3] 2 1 0
  //  !! Be careful, simplex are sorted by filtration value on insertion !!
  std::cout << "simplex_tree_insertion - first - 0" << std::endl;
  test_simplex_tree_contains(st, firstSimplex, 0); // (0) -> 0
  std::cout << "simplex_tree_insertion - second - 1" << std::endl;
  test_simplex_tree_contains(st, secondSimplex, 1); // (1) -> 1
  std::cout << "simplex_tree_insertion - third - 4" << std::endl;
  test_simplex_tree_contains(st, thirdSimplex, 4); // (0,1) -> 4
  std::cout << "simplex_tree_insertion - fourth - 2" << std::endl;
  test_simplex_tree_contains(st, fourthSimplex, 2); // (2) -> 2
  std::cout << "simplex_tree_insertion - fifth - 5" << std::endl;
  test_simplex_tree_contains(st, fifthSimplex, 5); // (2,0) -> 5
  std::cout << "simplex_tree_insertion - sixth - 6" << std::endl;
  test_simplex_tree_contains(st, sixthSimplex, 6); //(2,1) -> 6
  std::cout << "simplex_tree_insertion - seventh - 8" << std::endl;
  test_simplex_tree_contains(st, seventhSimplex, 8); // (2,1,0) -> 8
  std::cout << "simplex_tree_insertion - eighth - 3" << std::endl;
  test_simplex_tree_contains(st, eighthSimplex, 3); // (3) -> 3
  std::cout << "simplex_tree_insertion - nineth - 7" << std::endl;
  test_simplex_tree_contains(st, ninethSimplex, 7); // (3,0) -> 7

  // Display the Simplex_tree - Can not be done in the middle of 2 inserts
  std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
  std::cout << "   - dimension " << st.dimension() << std::endl;
  std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::cout << (int) vertex << " ";
    }
    std::cout << std::endl;
  }

}

BOOST_AUTO_TEST_CASE_TEMPLATE(NSimplexAndSubfaces_tree_insertion, typeST, list_of_tested_variants) {
  typedef std::pair<typename typeST::Simplex_handle, bool> typePairSimplexBool;
  typedef std::vector<typename typeST::Vertex_handle> typeVectorVertex;
  typedef std::pair<typeVectorVertex, typename typeST::Filtration_value> typeSimplex;
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF RECURSIVE INSERTION" << std::endl;
  typeST st;
  typePairSimplexBool returnValue;
  int position = 0;

  // ++ FIRST
  std::cout << "   - INSERT (2,1,0)" << std::endl;
  typeVectorVertex SimplexVector1{2, 1, 0};
  BOOST_CHECK(SimplexVector1.size() == 3);
  returnValue = st.insert_simplex_and_subfaces(SimplexVector1);

  BOOST_CHECK(st.num_vertices() == (size_t) 3); // +3 (2, 1 and 0 are not existing)

  // Check it is well inserted
  BOOST_CHECK(true == returnValue.second);
  position = 0;
  std::sort(SimplexVector1.begin(), SimplexVector1.end(), std::greater<typename typeST::Vertex_handle>());
  for (auto vertex : st.simplex_vertex_range(returnValue.first)) {
    // Check returned Simplex_handle
    std::cout << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector1[position] << std::endl;
    BOOST_CHECK(vertex == SimplexVector1[position]);
    position++;
  }

  // ++ SECOND
  std::cout << "   - INSERT 3" << std::endl;
  typeVectorVertex SimplexVector2{3};
  BOOST_CHECK(SimplexVector2.size() == 1);
  returnValue = st.insert_simplex_and_subfaces(SimplexVector2);

  BOOST_CHECK(st.num_vertices() == (size_t) 4); // +1 (3 is not existing)

  // Check it is well inserted
  BOOST_CHECK(true == returnValue.second);
  position = 0;
  std::sort(SimplexVector2.begin(), SimplexVector2.end(), std::greater<typename typeST::Vertex_handle>());
  for (auto vertex : st.simplex_vertex_range(returnValue.first)) {
    // Check returned Simplex_handle
    std::cout << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector2[position] << std::endl;
    BOOST_CHECK(vertex == SimplexVector2[position]);
    position++;
  }

  // ++ THIRD
  std::cout << "   - INSERT (0,3)" << std::endl;
  typeVectorVertex SimplexVector3{3, 0};
  BOOST_CHECK(SimplexVector3.size() == 2);
  returnValue = st.insert_simplex_and_subfaces(SimplexVector3);

  BOOST_CHECK(st.num_vertices() == (size_t) 4); // Not incremented (all are existing)

  // Check it is well inserted
  BOOST_CHECK(true == returnValue.second);
  position = 0;
  std::sort(SimplexVector3.begin(), SimplexVector3.end(), std::greater<typename typeST::Vertex_handle>());
  for (auto vertex : st.simplex_vertex_range(returnValue.first)) {
    // Check returned Simplex_handle
    std::cout << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector3[position] << std::endl;
    BOOST_CHECK(vertex == SimplexVector3[position]);
    position++;
  }

  // ++ FOURTH
  std::cout << "   - INSERT (1,0) (already inserted)" << std::endl;
  typeVectorVertex SimplexVector4{1, 0};
  BOOST_CHECK(SimplexVector4.size() == 2);
  returnValue = st.insert_simplex_and_subfaces(SimplexVector4);

  BOOST_CHECK(st.num_vertices() == (size_t) 4); // Not incremented (all are existing)

  // Check it was not inserted (already there from {2,1,0} insertion)
  BOOST_CHECK(false == returnValue.second);

  // ++ FIFTH
  std::cout << "   - INSERT (3,4,5)" << std::endl;
  typeVectorVertex SimplexVector5{3, 4, 5};
  BOOST_CHECK(SimplexVector5.size() == 3);
  returnValue = st.insert_simplex_and_subfaces(SimplexVector5);

  BOOST_CHECK(st.num_vertices() == (size_t) 6);

  // Check it is well inserted
  BOOST_CHECK(true == returnValue.second);
  position = 0;
  std::sort(SimplexVector5.begin(), SimplexVector5.end(), std::greater<typename typeST::Vertex_handle>());
  for (auto vertex : st.simplex_vertex_range(returnValue.first)) {
    // Check returned Simplex_handle
    std::cout << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector5[position] << std::endl;
    BOOST_CHECK(vertex == SimplexVector5[position]);
    position++;
  }

  // ++ SIXTH
  std::cout << "   - INSERT (0,1,6,7)" << std::endl;
  typeVectorVertex SimplexVector6{0, 1, 6, 7};
  BOOST_CHECK(SimplexVector6.size() == 4);
  returnValue = st.insert_simplex_and_subfaces(SimplexVector6);

  BOOST_CHECK(st.num_vertices() == (size_t) 8); // +2 (6 and 7 are not existing - 0 and 1 are already existing)

  // Check it is well inserted
  BOOST_CHECK(true == returnValue.second);
  position = 0;
  std::sort(SimplexVector6.begin(), SimplexVector6.end(), std::greater<typename typeST::Vertex_handle>());
  for (auto vertex : st.simplex_vertex_range(returnValue.first)) {
    // Check returned Simplex_handle
    std::cout << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector6[position] << std::endl;
    BOOST_CHECK(vertex == SimplexVector6[position]);
    position++;
  }
  
  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o           */
  /*  2   0   3\X/4           */
  /*            o             */
  /*            5             */
  /*                          */
  /* In other words:          */
  /*   A facet [2,1,0]        */
  /*   An edge [0,3]          */
  /*   A facet [3,4,5]        */
  /*   A cell  [0,1,6,7]      */

  typeSimplex simplexPair1 = std::make_pair(SimplexVector1, 0.0);
  typeSimplex simplexPair2 = std::make_pair(SimplexVector2, 0.0);
  typeSimplex simplexPair3 = std::make_pair(SimplexVector3, 0.0);
  typeSimplex simplexPair4 = std::make_pair(SimplexVector4, 0.0);
  typeSimplex simplexPair5 = std::make_pair(SimplexVector5, 0.0);
  typeSimplex simplexPair6 = std::make_pair(SimplexVector6, 0.0);
  test_simplex_tree_contains(st, simplexPair1, 6); // (2,1,0) is in position 6
  test_simplex_tree_contains(st, simplexPair2, 7); // (3) is in position 7
  test_simplex_tree_contains(st, simplexPair3, 8); // (3,0) is in position 8
  test_simplex_tree_contains(st, simplexPair4, 2); // (1,0) is in position 2
  test_simplex_tree_contains(st, simplexPair5, 14); // (3,4,5) is in position 14
  test_simplex_tree_contains(st, simplexPair6, 26); // (7,6,1,0) is in position 26

  // ------------------------------------------------------------------------------------------------------------------
  // Find in the simplex_tree
  // ------------------------------------------------------------------------------------------------------------------
  typeVectorVertex simpleSimplexVector{1};
  typename typeST::Simplex_handle simplexFound = st.find(simpleSimplexVector);
  std::cout << "**************IS THE SIMPLEX {1} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";
  // Check it is found
  BOOST_CHECK(simplexFound != st.null_simplex());

  typeVectorVertex unknownSimplexVector{15};
  simplexFound = st.find(unknownSimplexVector);
  std::cout << "**************IS THE SIMPLEX {15} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";
  // Check it is NOT found
  BOOST_CHECK(simplexFound == st.null_simplex());

  simplexFound = st.find(SimplexVector6);
  std::cout << "**************IS THE SIMPLEX {0,1,6,7} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";
  // Check it is found
  BOOST_CHECK(simplexFound != st.null_simplex());

  typeVectorVertex otherSimplexVector{1, 15};
  simplexFound = st.find(otherSimplexVector);
  std::cout << "**************IS THE SIMPLEX {15,1} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";
  // Check it is NOT found
  BOOST_CHECK(simplexFound == st.null_simplex());

  typeVectorVertex invSimplexVector{1, 2, 0};
  simplexFound = st.find(invSimplexVector);
  std::cout << "**************IS THE SIMPLEX {1,2,0} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::cout << "***+ YES IT IS!\n";
  else
    std::cout << "***- NO IT ISN'T\n";
  // Check it is found
  BOOST_CHECK(simplexFound != st.null_simplex());

  // Display the Simplex_tree - Can not be done in the middle of 2 inserts
  std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
  std::cout << "   - dimension " << st.dimension() << std::endl;
  std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::cout << (int) vertex << " ";
    }
    std::cout << std::endl;
  }
}

template<class typeST, class Vertex_handle>
void test_cofaces(typeST& st, const std::vector<Vertex_handle>& expected, int dim, const std::vector<typename typeST::Simplex_handle>& res) {
  typename typeST::Cofaces_simplex_range cofaces;
  if (dim == 0)
    cofaces = st.star_simplex_range(st.find(expected));
  else
    cofaces = st.cofaces_simplex_range(st.find(expected), dim);
  for (auto simplex = cofaces.begin(); simplex != cofaces.end(); ++simplex) {
    typename typeST::Simplex_vertex_range rg = st.simplex_vertex_range(*simplex);
    for (auto vertex = rg.begin(); vertex != rg.end(); ++vertex) {
      std::cout << "(" << *vertex << ")";
    }
    std::cout << std::endl;
    BOOST_CHECK(std::find(res.begin(), res.end(), *simplex) != res.end());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(coface_on_simplex_tree, typeST, list_of_tested_variants) {
  typedef std::vector<typename typeST::Vertex_handle> typeVectorVertex;
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST COFACE ALGORITHM" << std::endl;
  typeST st;

  typeVectorVertex SimplexVector{2, 1, 0};
  st.insert_simplex_and_subfaces(SimplexVector);

  SimplexVector = {3, 0};
  st.insert_simplex_and_subfaces(SimplexVector);

  SimplexVector = {3, 4, 5};
  st.insert_simplex_and_subfaces(SimplexVector);

  SimplexVector = {0, 1, 6, 7};
  st.insert_simplex_and_subfaces(SimplexVector);

  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o           */
  /*  2   0   3\X/4           */
  /*            o             */
  /*            5             */

  std::vector<typename typeST::Vertex_handle> simplex_result;
  std::vector<typename typeST::Simplex_handle> result;
  std::cout << "First test - Star of (3):" << std::endl;

  simplex_result = {3};
  result.push_back(st.find(simplex_result));

  simplex_result = {3, 0};
  result.push_back(st.find(simplex_result));

  simplex_result = {4, 3};
  result.push_back(st.find(simplex_result));

  simplex_result = {5, 4, 3};
  result.push_back(st.find(simplex_result));

  simplex_result = {5, 3};
  result.push_back(st.find(simplex_result));
  simplex_result.clear();

  std::vector<typename typeST::Vertex_handle> vertex = {3};
  test_cofaces(st, vertex, 0, result);
  vertex.clear();
  result.clear();

  vertex.push_back(1);
  vertex.push_back(7);
  std::cout << "Second test - Star of (1,7): " << std::endl;

  simplex_result = {7, 1};
  result.push_back(st.find(simplex_result));

  simplex_result = {7, 6, 1, 0};
  result.push_back(st.find(simplex_result));

  simplex_result = {7, 1, 0};
  result.push_back(st.find(simplex_result));

  simplex_result = {7, 6, 1};
  result.push_back(st.find(simplex_result));

  test_cofaces(st, vertex, 0, result);
  result.clear();

  std::cout << "Third test - 2-dimension Cofaces of simplex(1,7) : " << std::endl;

  simplex_result = {7, 1, 0};
  result.push_back(st.find(simplex_result));

  simplex_result = {7, 6, 1};
  result.push_back(st.find(simplex_result));

  test_cofaces(st, vertex, 1, result);
  result.clear();

  std::cout << "Cofaces with a codimension too high (codimension + vetices > tree.dimension) :" << std::endl;
  test_cofaces(st, vertex, 5, result);

  //std::cout << "Cofaces with an empty codimension" << std::endl;
  //test_cofaces(st, vertex, -1, result);
  //    std::cout << "Cofaces in an empty simplex tree" << std::endl;
  //   typeST empty_tree;
  //    test_cofaces(empty_tree, vertex, 1, result);
  //std::cout << "Cofaces of an empty simplex" << std::endl;
  //vertex.clear();
  // test_cofaces(st, vertex, 1, result);

}

BOOST_AUTO_TEST_CASE_TEMPLATE(copy_move_on_simplex_tree, typeST, list_of_tested_variants) {
  typedef std::vector<typename typeST::Vertex_handle> typeVectorVertex;
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST COPY MOVE CONSTRUCTORS" << std::endl;
  typeST st;

  typeVectorVertex SimplexVector{2, 1, 0};
  st.insert_simplex_and_subfaces(SimplexVector);

  SimplexVector = {3, 0};
  st.insert_simplex_and_subfaces(SimplexVector);

  SimplexVector = {3, 4, 5};
  st.insert_simplex_and_subfaces(SimplexVector);

  SimplexVector = {0, 1, 6, 7};
  st.insert_simplex_and_subfaces(SimplexVector);

  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o           */
  /*  2   0   3\X/4           */
  /*            o             */
  /*            5             */

  std::cout << "Printing st - address = " << &st << std::endl;

  // Copy constructor  
  typeST st_copy = st;
  std::cout << "Printing a copy of st - address = " << &st_copy << std::endl;

  // Check the data are the same
  BOOST_CHECK(st == st_copy);
  // Check there is a new simplex tree reference
  BOOST_CHECK(&st != &st_copy);

  // Move constructor  
  typeST st_move = std::move(st);
  std::cout << "Printing a move of st - address = " << &st_move << std::endl;

  // Check the data are the same
  BOOST_CHECK(st_move == st_copy);
  // Check there is a new simplex tree reference
  BOOST_CHECK(&st_move != &st_copy);
  BOOST_CHECK(&st_move != &st);
  
  typeST st_empty;
  // Check st has been emptied by the move
  BOOST_CHECK(st == st_empty);
  BOOST_CHECK(st.dimension() == -1);
  BOOST_CHECK(st.num_simplices() == 0);
  BOOST_CHECK(st.num_vertices() == (size_t)0);
  
  std::cout << "Printing st once again- address = " << &st << std::endl;
}

template<class typeST>
void test_simplex_is_vertex(typeST& st, typename typeST::Simplex_handle sh, typename typeST::Vertex_handle v) {
  BOOST_CHECK(st.dimension(sh) == 0);
  auto&& r = st.simplex_vertex_range(sh);
  auto i = std::begin(r);
  BOOST_CHECK(*i == v);
  BOOST_CHECK(++i == std::end(r));
}

BOOST_AUTO_TEST_CASE(non_contiguous) {
  typedef Simplex_tree<> typeST;
  typedef typeST::Simplex_handle Simplex_handle;
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST NON-CONTIGUOUS VERTICES" << std::endl;
  typeST st;
  typeST::Vertex_handle e[] = {3,-7};
  std::cout << "Insert" << std::endl;
  st.insert_simplex_and_subfaces(e);
  BOOST_CHECK(st.num_vertices() == 2);
  BOOST_CHECK(st.num_simplices() == 3);
  std::cout << "Find" << std::endl;
  Simplex_handle sh = st.find(e);
  BOOST_CHECK(sh != st.null_simplex());
  std::cout << "Endpoints" << std::endl;
  auto p = st.endpoints(sh);
  test_simplex_is_vertex(st, p.first, 3);
  test_simplex_is_vertex(st, p.second, -7);
  std::cout << "Boundary" << std::endl;
  auto&& b = st.boundary_simplex_range(sh);
  auto i = std::begin(b);
  test_simplex_is_vertex(st, *i, -7);
  test_simplex_is_vertex(st, *++i, 3);
  BOOST_CHECK(++i == std::end(b));
}

BOOST_AUTO_TEST_CASE(make_filtration_non_decreasing) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "MAKE FILTRATION NON DECREASING" << std::endl;
  typedef Simplex_tree<> typeST;
  typeST st;

  st.insert_simplex_and_subfaces({2, 1, 0}, 2.0);
  st.insert_simplex_and_subfaces({3, 0}, 2.0);
  st.insert_simplex_and_subfaces({3, 4, 5}, 2.0);
  
  /* Inserted simplex:     */
  /*    1                  */
  /*    o                  */
  /*   /X\                 */
  /*  o---o---o---o        */
  /*  2   0   3\X/4        */
  /*            o          */
  /*            5          */

  std::cout << "Check default insertion ensures the filtration values are non decreasing" << std::endl;
  BOOST_CHECK(!st.make_filtration_non_decreasing());

  // Because of non decreasing property of simplex tree, { 0 } , { 1 } and { 0, 1 } are going to be set from value 2.0
  // to 1.0
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, 1.0);
  
  // Inserted simplex:
  //    1   6
  //    o---o
  //   /X\7/
  //  o---o---o---o
  //  2   0   3\X/4
  //            o
  //            5
  
  std::cout << "Check default second insertion ensures the filtration values are non decreasing" << std::endl;
  BOOST_CHECK(!st.make_filtration_non_decreasing());
  
  // Copy original simplex tree
  typeST st_copy = st;

  // Modify specific values for st to become like st_copy thanks to make_filtration_non_decreasing
  st.assign_filtration(st.find({0,1,6,7}), 0.8);
  st.assign_filtration(st.find({0,1,6}), 0.9);
  st.assign_filtration(st.find({0,6}), 0.6);
  st.assign_filtration(st.find({3,4,5}), 1.2);
  st.assign_filtration(st.find({3,4}), 1.1);
  st.assign_filtration(st.find({4,5}), 1.99);
  
  std::cout << "Check the simplex_tree is rolled back in case of decreasing filtration values" << std::endl;
  BOOST_CHECK(st.make_filtration_non_decreasing());
  BOOST_CHECK(st == st_copy);

  // Other simplex tree
  typeST st_other;
  st_other.insert_simplex_and_subfaces({2, 1, 0}, 3.0);  // This one is different from st
  st_other.insert_simplex_and_subfaces({3, 0}, 2.0);
  st_other.insert_simplex_and_subfaces({3, 4, 5}, 2.0);
  st_other.insert_simplex_and_subfaces({0, 1, 6, 7}, 1.0);

  // Modify specific values for st to become like st_other thanks to make_filtration_non_decreasing
  st.assign_filtration(st.find({2}), 3.0);
  // By modifying just the simplex {2}
  // {0,1,2}, {1,2} and {0,2} will be modified
  
  std::cout << "Check the simplex_tree is repaired in case of decreasing filtration values" << std::endl;
  BOOST_CHECK(st.make_filtration_non_decreasing());
  BOOST_CHECK(st == st_other);

  // Modify specific values for st still to be non-decreasing
  st.assign_filtration(st.find({0,1,2}), 10.0);
  st.assign_filtration(st.find({0,2}), 9.0);
  st.assign_filtration(st.find({0,1,6,7}), 50.0);
  st.assign_filtration(st.find({0,1,6}), 49.0);
  st.assign_filtration(st.find({0,1,7}), 48.0);
  // Other copy simplex tree
  typeST st_other_copy = st;
  
  std::cout << "Check the simplex_tree is not modified in case of non-decreasing filtration values" << std::endl;
  BOOST_CHECK(!st.make_filtration_non_decreasing());
  BOOST_CHECK(st == st_other_copy);
  
}


typedef boost::mpl::list<boost::adjacency_list<boost::setS, boost::vecS, boost::directedS,
                                               boost::property<vertex_filtration_t, double>,
                                               boost::property<edge_filtration_t, double>>,
                         boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                               boost::property<vertex_filtration_t, double>,
                                               boost::property<edge_filtration_t, double>>,
                         boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS,
                                               boost::property<vertex_filtration_t, double>,
                                               boost::property<edge_filtration_t, double>>,
                         boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                                               boost::property<vertex_filtration_t, double>,
                                               boost::property<edge_filtration_t, double>>,
                         boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                               boost::property<vertex_filtration_t, double>,
                                               boost::property<edge_filtration_t, double>>,
                         boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
                                               boost::property<vertex_filtration_t, double>,
                                               boost::property<edge_filtration_t, double>>> list_of_graph_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_insert_graph, Graph, list_of_graph_variants) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "INSERT GRAPH" << std::endl;

  Graph g(3);
  // filtration value 0 everywhere
  put(Gudhi::vertex_filtration_t(), g, 0, 0);
  put(Gudhi::vertex_filtration_t(), g, 1, 0);
  put(Gudhi::vertex_filtration_t(), g, 2, 0);
  // vertices don't always occur in sorted order
  add_edge(0, 1, 1.1, g);
  add_edge(2, 0, 2.2, g);
  add_edge(2, 1, 3.3, g);

  Simplex_tree<> st1;
  st1.insert_graph(g);
  BOOST_CHECK(st1.num_simplices() == 6);

  // edges can have multiplicity in the graph unless we replace the first vecS with (hash_)setS
  add_edge(1, 0, 1.1, g);
  add_edge(1, 2, 3.3, g);
  add_edge(0, 2, 2.2, g);
  add_edge(0, 1, 1.1, g);
  add_edge(2, 1, 3.3, g);
  add_edge(2, 0, 2.2, g);
  Simplex_tree<> st2;
  st2.insert_graph(g);
  BOOST_CHECK(st2.num_simplices() == 6);

  std::cout << "st1 is" << std::endl;
  std::cout << st1 << std::endl;

  std::cout << "st2 is" << std::endl;
  std::cout << st2 << std::endl;

  BOOST_CHECK(st1 == st2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(insert_duplicated_vertices, typeST, list_of_tested_variants) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST INSERT DUPLICATED VERTICES" << std::endl;
  typeST st;

  typename typeST::Simplex_handle sh;
  bool success = false;
  std::tie(sh, success) = st.insert_simplex_and_subfaces({1});
  BOOST_CHECK(success);
  BOOST_CHECK(sh != st.null_simplex());
  std::cout << "st.dimension(sh)= " << st.dimension(sh) << std::endl;
  BOOST_CHECK(st.dimension(sh) == 0);
  std::tie(sh, success) = st.insert_simplex_and_subfaces({2, 2});
  BOOST_CHECK(success);
  BOOST_CHECK(sh != st.null_simplex());
  std::cout << "st.dimension(sh)= " << st.dimension(sh) << std::endl;
  BOOST_CHECK(st.dimension(sh) == 0);
  std::tie(sh, success) = st.insert_simplex_and_subfaces({3, 3, 3});
  BOOST_CHECK(success);
  BOOST_CHECK(sh != st.null_simplex());
  std::cout << "st.dimension(sh)= " << st.dimension(sh) << std::endl;
  BOOST_CHECK(st.dimension(sh) == 0);
  std::tie(sh, success) = st.insert_simplex_and_subfaces({4, 4, 4, 4});
  BOOST_CHECK(success);
  BOOST_CHECK(sh != st.null_simplex());
  std::cout << "st.dimension(sh)= " << st.dimension(sh) << std::endl;
  BOOST_CHECK(st.dimension(sh) == 0);

  std::cout << "dimension =" << st.dimension() << " - num_vertices = " << st.num_vertices()
            << " - num_simplices = " << st.num_simplices() << std::endl;
  BOOST_CHECK(st.dimension() == 0);
  BOOST_CHECK(st.num_simplices() == st.num_vertices());

  std::tie(sh, success) = st.insert_simplex_and_subfaces({2, 1, 1, 2});
  BOOST_CHECK(success);
  BOOST_CHECK(sh != st.null_simplex());
  std::cout << "st.dimension(sh)= " << st.dimension(sh) << std::endl;
  BOOST_CHECK(st.dimension(sh) == 1);

  std::cout << "dimension =" << st.dimension() << " - num_vertices = " << st.num_vertices()
            << " - num_simplices = " << st.num_simplices() << std::endl;
  BOOST_CHECK(st.dimension() == 1);
  BOOST_CHECK(st.num_simplices() == st.num_vertices() + 1);

  // Already inserted
  std::tie(sh, success) = st.insert_simplex_and_subfaces({1, 2, 2, 1});
  BOOST_CHECK(!success);
  BOOST_CHECK(sh == st.null_simplex());

  std::cout << "dimension =" << st.dimension() << " - num_vertices = " << st.num_vertices()
            << " - num_simplices = " << st.num_simplices() << std::endl;
  BOOST_CHECK(st.dimension() == 1);
  BOOST_CHECK(st.num_simplices() == st.num_vertices() + 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(generators, typeST, list_of_tested_variants) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST FIND GENERATORS" << std::endl;
  {
    typeST st;
    st.insert_simplex_and_subfaces({0,1,2,3,4,5,6},0);
    st.assign_filtration(st.find({0,2,4}), 10);
    st.assign_filtration(st.find({1,5}), 20);
    st.assign_filtration(st.find({1,2,4}), 30);
    st.assign_filtration(st.find({3}), 5);
    st.make_filtration_non_decreasing();
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
    typeST st;
    st.insert_simplex_and_subfaces({0,1}, 8);
    st.insert_simplex_and_subfaces({0,2}, 10);
    st.insert_simplex_and_subfaces({3,4}, 6);
    st.insert_simplex_and_subfaces({1,2}, 5);
    st.insert_simplex_and_subfaces({1,5}, 4);
    st.insert_simplex_and_subfaces({0,5}, 3);
    st.insert_simplex_and_subfaces({2,5}, 2);
    st.insert_simplex_and_subfaces({1,3}, 9);
    st.expansion(50);
    BOOST_CHECK(st.edge_with_same_filtration(st.find({0,1,2,5}))==st.find({0,2}));
    BOOST_CHECK(st.edge_with_same_filtration(st.find({1,5}))==st.find({1,5}));
  }
}
