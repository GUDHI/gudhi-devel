/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux, Vincent Rouvreau
 *
 *    Copyright (C) 2023 Inria
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
#include <iterator>  // for std::distance
#include <cstddef>  // for std::size_t

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_multi"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree_multi.h>
#include <gudhi/One_critical_filtration.h>
#include <gudhi/Multi_critical_filtration.h>

using namespace Gudhi;
using OneCriticalFiltration = Gudhi::multi_filtration::One_critical_filtration<float>;
using Stree_options_multi_filtration = Gudhi::multi_persistence::Simplex_tree_options_multidimensional_filtration<OneCriticalFiltration>;

using typeST_STD = Simplex_tree<Simplex_tree_options_full_featured>;
using Stree_multi = Simplex_tree<Stree_options_multi_filtration>;

typedef boost::mpl::list<Stree_multi> list_of_tested_variants;

template<class typeST>
void test_empty_simplex_tree(typeST& tst) {
  typedef typename typeST::Vertex_handle Vertex_handle;
  const Vertex_handle DEFAULT_VERTEX_VALUE = Vertex_handle(- 1);
  BOOST_CHECK(tst.null_vertex() == DEFAULT_VERTEX_VALUE);
  BOOST_CHECK(tst.num_vertices() == (size_t) 0);
  BOOST_CHECK(tst.num_simplices() == (size_t) 0);
  BOOST_CHECK(tst.is_empty());
  BOOST_CHECK(tst.num_simplices_by_dimension() == std::vector<size_t>());
  typename typeST::Siblings* STRoot = tst.root();
  BOOST_CHECK(STRoot != nullptr);
  BOOST_CHECK(STRoot->oncles() == nullptr);
  BOOST_CHECK(STRoot->parent() == DEFAULT_VERTEX_VALUE);
  BOOST_CHECK(tst.dimension() == -1);
}

template<class typeST>
void test_iterators_on_empty_simplex_tree(typeST& tst) {
  std::clog << "Iterator on vertices: " << std::endl;
  for (auto vertex : tst.complex_vertex_range()) {
    std::clog << "vertice:" << vertex << std::endl;
    BOOST_CHECK(false); // shall be empty
  }
  std::clog << "Iterator on simplices: " << std::endl;
  for (auto simplex : tst.complex_simplex_range()) {
    BOOST_CHECK(simplex != simplex); // shall be empty - to remove warning of non-used simplex
  }

  std::clog
      << "Iterator on Simplices in the filtration, with [filtration value]:"
      << std::endl;
  for (auto f_simplex : tst.filtration_simplex_range()) {
    BOOST_CHECK(false); // shall be empty
    std::clog << "test_iterators_on_empty_simplex_tree - filtration="
        << tst.filtration(f_simplex) << std::endl;
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_when_empty, typeST, list_of_tested_variants) {
  typedef std::pair<typename typeST::Simplex_handle, bool> typePairSimplexBool;
  typedef std::vector<typename typeST::Vertex_handle> typeVectorVertex;

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF DEFAULT CONSTRUCTOR" << std::endl;
  typeST st;

  test_empty_simplex_tree(st);

  test_iterators_on_empty_simplex_tree(st);
  // TEST OF EMPTY INSERTION
  std::clog << "TEST OF EMPTY INSERTION" << std::endl;
  typeVectorVertex simplexVectorEmpty;
  BOOST_CHECK(simplexVectorEmpty.empty() == true);
  typePairSimplexBool returnEmptyValue = st.insert_simplex(simplexVectorEmpty, 0.0);
  BOOST_CHECK(returnEmptyValue.first == typeST::null_simplex());
  BOOST_CHECK(returnEmptyValue.second == true);

  test_empty_simplex_tree(st);

  test_iterators_on_empty_simplex_tree(st);
}

bool AreAlmostTheSame(Gudhi::multi_filtration::One_critical_filtration<float> a,
                      Gudhi::multi_filtration::One_critical_filtration<float> b) {
  assert(a.size() == b.size());
  for (auto i=0u; i<a.size();i++)
    if (std::fabs(a[i] - b[i]) > std::numeric_limits<float>::epsilon())
      return false;
  return true;
}


template<class typeST, class typeSimplex>
void test_simplex_tree_contains(typeST& simplexTree, typeSimplex& simplex, int pos) {
  auto f_simplex = simplexTree.filtration_simplex_range().begin() + pos;

  std::clog << "test_simplex_tree_contains - filtration=" << simplexTree.filtration(*f_simplex) << "||" << simplex.second << std::endl;
  BOOST_CHECK(AreAlmostTheSame(simplexTree.filtration(*f_simplex), simplex.second));

  int simplexIndex = simplex.first.size() - 1;
  std::sort(simplex.first.begin(), simplex.first.end()); // if the simplex wasn't sorted, the next test could fail
  for (auto vertex : simplexTree.simplex_vertex_range(*f_simplex)) {
    std::clog << "test_simplex_tree_contains - vertex=" << vertex << "||" << simplex.first.at(simplexIndex) << std::endl;
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
  BOOST_CHECK(shReturned != typeST::null_simplex());
}

// Global variables
int dim_max = -1;

template<class typeST, class Filtration_value>
void set_and_test_simplex_tree_dim_fil(typeST& simplexTree, int vectorSize, const Filtration_value& fil) {
  if (vectorSize > dim_max + 1) {
    dim_max = vectorSize - 1;
    simplexTree.set_dimension(dim_max);
    std::clog << "   set_and_test_simplex_tree_dim_fil - dim_max=" << dim_max
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
  const Filtration_value FIRST_FILTRATION_VALUE = {0.1,0.15};
  const Filtration_value SECOND_FILTRATION_VALUE = {0.2, 0.25};
  const Filtration_value THIRD_FILTRATION_VALUE = {0.3, 0.35};
  const Filtration_value FOURTH_FILTRATION_VALUE = {0.4, 0.45};
  // reset since we run the test several times
  dim_max = -1;

  // TEST OF INSERTION
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF INSERTION" << std::endl;
  typeST st;

  // ++ FIRST
  std::clog << "   - INSERT 0" << std::endl;
  typeVectorVertex firstSimplexVector{0};
  BOOST_CHECK(firstSimplexVector.size() == 1);
  typeSimplex firstSimplex = std::make_pair(firstSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
  typePairSimplexBool returnValue = st.insert_simplex(firstSimplex.first, firstSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, firstSimplexVector.size(), firstSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 1);

  // ++ SECOND
  std::clog << "   - INSERT 1" << std::endl;
  typeVectorVertex secondSimplexVector{1};
  BOOST_CHECK(secondSimplexVector.size() == 1);
  typeSimplex secondSimplex = std::make_pair(secondSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
  returnValue = st.insert_simplex(secondSimplex.first, secondSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, secondSimplexVector.size(), secondSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 2);

  // ++ THIRD
  std::clog << "   - INSERT (0,1)" << std::endl;
  typeVectorVertex thirdSimplexVector{0, 1};
  BOOST_CHECK(thirdSimplexVector.size() == 2);
  typeSimplex thirdSimplex = std::make_pair(thirdSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
  returnValue = st.insert_simplex(thirdSimplex.first, thirdSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, thirdSimplexVector.size(), thirdSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 2); // Not incremented !!

  // ++ FOURTH
  std::clog << "   - INSERT 2" << std::endl;
  typeVectorVertex fourthSimplexVector{2};
  BOOST_CHECK(fourthSimplexVector.size() == 1);
  typeSimplex fourthSimplex = std::make_pair(fourthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
  returnValue = st.insert_simplex(fourthSimplex.first, fourthSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, fourthSimplexVector.size(), fourthSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 3);

  // ++ FIFTH
  std::clog << "   - INSERT (2,0)" << std::endl;
  typeVectorVertex fifthSimplexVector{2, 0};
  BOOST_CHECK(fifthSimplexVector.size() == 2);
  typeSimplex fifthSimplex = std::make_pair(fifthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
  returnValue = st.insert_simplex(fifthSimplex.first, fifthSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, fifthSimplexVector.size(), fifthSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 3); // Not incremented !!

  // ++ SIXTH
  std::clog << "   - INSERT (2,1)" << std::endl;
  typeVectorVertex sixthSimplexVector{2, 1};
  BOOST_CHECK(sixthSimplexVector.size() == 2);
  typeSimplex sixthSimplex = std::make_pair(sixthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
  returnValue = st.insert_simplex(sixthSimplex.first, sixthSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, sixthSimplexVector.size(), sixthSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 3); // Not incremented !!

  // ++ SEVENTH
  std::clog << "   - INSERT (2,1,0)" << std::endl;
  typeVectorVertex seventhSimplexVector{2, 1, 0};
  BOOST_CHECK(seventhSimplexVector.size() == 3);
  typeSimplex seventhSimplex = std::make_pair(seventhSimplexVector, Filtration_value(THIRD_FILTRATION_VALUE));
  returnValue = st.insert_simplex(seventhSimplex.first, seventhSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, seventhSimplexVector.size(), seventhSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 3); // Not incremented !!

  // ++ EIGHTH
  std::clog << "   - INSERT 3" << std::endl;
  typeVectorVertex eighthSimplexVector{3};
  BOOST_CHECK(eighthSimplexVector.size() == 1);
  typeSimplex eighthSimplex = std::make_pair(eighthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
  returnValue = st.insert_simplex(eighthSimplex.first, eighthSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, eighthSimplexVector.size(), eighthSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 4);

  // ++ NINTH
  std::clog << "   - INSERT (3,0)" << std::endl;
  typeVectorVertex ninethSimplexVector{3, 0};
  BOOST_CHECK(ninethSimplexVector.size() == 2);
  typeSimplex ninethSimplex = std::make_pair(ninethSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
  returnValue = st.insert_simplex(ninethSimplex.first, ninethSimplex.second);

  test_simplex_tree_insert_returns_true<typeST>(returnValue);
  set_and_test_simplex_tree_dim_fil(st, ninethSimplexVector.size(), ninethSimplex.second);
  BOOST_CHECK(st.num_vertices() == (size_t) 4); // Not incremented !!

  // ++ TENTH
  std::clog << "   - INSERT 0 (already inserted)" << std::endl;
  typeVectorVertex tenthSimplexVector{0};
  BOOST_CHECK(tenthSimplexVector.size() == 1);
  // With a different filtration value
  typeSimplex tenthSimplex = std::make_pair(tenthSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE));
  returnValue = st.insert_simplex(tenthSimplex.first, tenthSimplex.second);

  BOOST_CHECK(returnValue.second == false);
  // Simplex_handle = boost::container::flat_map< typeST::Vertex_handle, Node >::iterator
  typename typeST::Simplex_handle shReturned = returnValue.first;
  BOOST_CHECK(shReturned == typeST::null_simplex());
  std::clog << "st.num_vertices()=" << st.num_vertices() << std::endl;
  BOOST_CHECK(st.num_vertices() == (size_t) 4); // Not incremented !!
  BOOST_CHECK(st.dimension() == dim_max);

  // ++ ELEVENTH
  std::clog << "   - INSERT (2,1,0) (already inserted)" << std::endl;
  typeVectorVertex eleventhSimplexVector{2, 1, 0};
  BOOST_CHECK(eleventhSimplexVector.size() == 3);
  typeSimplex eleventhSimplex = std::make_pair(eleventhSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE));
  returnValue = st.insert_simplex(eleventhSimplex.first, eleventhSimplex.second);

  BOOST_CHECK(returnValue.second == false);
  // Simplex_handle = boost::container::flat_map< typeST::Vertex_handle, Node >::iterator
  shReturned = returnValue.first;
  BOOST_CHECK(shReturned == typeST::null_simplex());
  BOOST_CHECK(st.num_vertices() == (size_t) 4); // Not incremented !!
  BOOST_CHECK(st.dimension() == dim_max);
  BOOST_CHECK(st.num_simplices_by_dimension() == std::vector<size_t>({4, 4, 1}));

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
  std::clog << "simplex_tree_insertion - first - 0" << std::endl;
  test_simplex_tree_contains(st, firstSimplex, 0); // (0) -> 0
  std::clog << "simplex_tree_insertion - second - 1" << std::endl;
  test_simplex_tree_contains(st, secondSimplex, 1); // (1) -> 1
  std::clog << "simplex_tree_insertion - third - 4" << std::endl;
  test_simplex_tree_contains(st, thirdSimplex, 4); // (0,1) -> 4
  std::clog << "simplex_tree_insertion - fourth - 2" << std::endl;
  test_simplex_tree_contains(st, fourthSimplex, 2); // (2) -> 2
  std::clog << "simplex_tree_insertion - fifth - 5" << std::endl;
  test_simplex_tree_contains(st, fifthSimplex, 5); // (2,0) -> 5
  std::clog << "simplex_tree_insertion - sixth - 6" << std::endl;
  test_simplex_tree_contains(st, sixthSimplex, 6); //(2,1) -> 6
  std::clog << "simplex_tree_insertion - seventh - 8" << std::endl;
  test_simplex_tree_contains(st, seventhSimplex, 8); // (2,1,0) -> 8
  std::clog << "simplex_tree_insertion - eighth - 3" << std::endl;
  test_simplex_tree_contains(st, eighthSimplex, 3); // (3) -> 3
  std::clog << "simplex_tree_insertion - ninth - 7" << std::endl;
  test_simplex_tree_contains(st, ninethSimplex, 7); // (3,0) -> 7

  // Display the Simplex_tree - Can not be done in the middle of 2 inserts
  std::clog << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
  std::clog << "   - dimension " << st.dimension() << std::endl;
  std::clog << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::clog << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::clog << (int) vertex << " ";
    }
    std::clog << std::endl;
  }

}

BOOST_AUTO_TEST_CASE_TEMPLATE(NSimplexAndSubfaces_tree_insertion, typeST, list_of_tested_variants) {
  typedef std::pair<typename typeST::Simplex_handle, bool> typePairSimplexBool;
  typedef std::vector<typename typeST::Vertex_handle> typeVectorVertex;
  typedef typename typeST::Filtration_value F;
  typedef std::pair<typeVectorVertex, F> typeSimplex;
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF RECURSIVE INSERTION" << std::endl;
  typeST st;
  typePairSimplexBool returnValue;
  int position = 0;

  // ++ FIRST
  std::clog << "   - INSERT (2,1,0)" << std::endl;
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
    std::clog << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector1[position] << std::endl;
    BOOST_CHECK(vertex == SimplexVector1[position]);
    position++;
  }

  // ++ SECOND
  std::clog << "   - INSERT 3" << std::endl;
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
    std::clog << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector2[position] << std::endl;
    BOOST_CHECK(vertex == SimplexVector2[position]);
    position++;
  }

  // ++ THIRD
  std::clog << "   - INSERT (0,3)" << std::endl;
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
    std::clog << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector3[position] << std::endl;
    BOOST_CHECK(vertex == SimplexVector3[position]);
    position++;
  }

  // ++ FOURTH
  std::clog << "   - INSERT (1,0) (already inserted)" << std::endl;
  typeVectorVertex SimplexVector4{1, 0};
  BOOST_CHECK(SimplexVector4.size() == 2);
  returnValue = st.insert_simplex_and_subfaces(SimplexVector4);

  BOOST_CHECK(st.num_vertices() == (size_t) 4); // Not incremented (all are existing)

  // Check it was not inserted (already there from {2,1,0} insertion)
  BOOST_CHECK(false == returnValue.second);

  // ++ FIFTH
  std::clog << "   - INSERT (3,4,5)" << std::endl;
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
    std::clog << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector5[position] << std::endl;
    BOOST_CHECK(vertex == SimplexVector5[position]);
    position++;
  }

  // ++ SIXTH
  std::clog << "   - INSERT (0,1,6,7)" << std::endl;
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
    std::clog << "vertex = " << vertex << " | vector[" << position << "] = " << SimplexVector6[position] << std::endl;
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

  typeSimplex simplexPair1 = std::make_pair(SimplexVector1, F());
  typeSimplex simplexPair2 = std::make_pair(SimplexVector2, F());
  typeSimplex simplexPair3 = std::make_pair(SimplexVector3, F());
  typeSimplex simplexPair4 = std::make_pair(SimplexVector4, F());
  typeSimplex simplexPair5 = std::make_pair(SimplexVector5, F());
  typeSimplex simplexPair6 = std::make_pair(SimplexVector6, F());
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
  std::clog << "**************IS THE SIMPLEX {1} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";
  // Check it is found
  BOOST_CHECK(simplexFound != st.null_simplex());

  typeVectorVertex unknownSimplexVector{15};
  simplexFound = st.find(unknownSimplexVector);
  std::clog << "**************IS THE SIMPLEX {15} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";
  // Check it is NOT found
  BOOST_CHECK(simplexFound == st.null_simplex());

  simplexFound = st.find(SimplexVector6);
  std::clog << "**************IS THE SIMPLEX {0,1,6,7} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";
  // Check it is found
  BOOST_CHECK(simplexFound != st.null_simplex());

  typeVectorVertex otherSimplexVector{1, 15};
  simplexFound = st.find(otherSimplexVector);
  std::clog << "**************IS THE SIMPLEX {15,1} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";
  // Check it is NOT found
  BOOST_CHECK(simplexFound == st.null_simplex());

  typeVectorVertex invSimplexVector{1, 2, 0};
  simplexFound = st.find(invSimplexVector);
  std::clog << "**************IS THE SIMPLEX {1,2,0} IN THE SIMPLEX TREE ?\n";
  if (simplexFound != st.null_simplex())
    std::clog << "***+ YES IT IS!\n";
  else
    std::clog << "***- NO IT ISN'T\n";
  // Check it is found
  BOOST_CHECK(simplexFound != st.null_simplex());

  // Display the Simplex_tree - Can not be done in the middle of 2 inserts
  std::clog << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
  std::clog << "   - dimension " << st.dimension() << std::endl;
  std::clog << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::clog << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::clog << (int) vertex << " ";
    }
    std::clog << std::endl;
  }
}


template<class typeST>
void test_simplex_is_vertex(typeST& st, typename typeST::Simplex_handle sh, typename typeST::Vertex_handle v) {
  BOOST_CHECK(st.dimension(sh) == 0);
  auto&& r = st.simplex_vertex_range(sh);
  auto i = std::begin(r);
  BOOST_CHECK(*i == v);
  BOOST_CHECK(++i == std::end(r));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_reset_filtration, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST RESET FILTRATION" << std::endl;
  typeST st;

  st.insert_simplex_and_subfaces({2, 1, 0}, OneCriticalFiltration({2.,1.}));
  st.insert_simplex_and_subfaces({3, 0}, OneCriticalFiltration({1.,2.}));
  st.insert_simplex_and_subfaces({3, 4, 5}, OneCriticalFiltration({3.,4.}));
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, OneCriticalFiltration({4.,3.}));
  std::cout <<"TRUC "<< st.filtration(st.find({2,1,0})) << std::endl;
  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o           */
  /*  2   0   3\X/4           */
  /*            o             */
  /*            5             */

  for (auto f_simplex : st.skeleton_simplex_range(3)) {
    std::clog << "vertex = (";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << ",";
    }
    std::clog << ") - filtration = " << st.filtration(f_simplex);
    std::clog << " - dimension = " << st.dimension(f_simplex) << std::endl;
    // Guaranteed by construction
    BOOST_CHECK(st.filtration(f_simplex) >= OneCriticalFiltration({1.,1.}));
  }

  // dimension until 5 even if simplex tree is of dimension 3 to test the limits
  for(int dimension = 5; dimension >= 0; dimension --) {
    std::clog << "### reset_filtration - dimension = " << dimension << "\n";
    st.reset_filtration(st.inf_, dimension);
    for (auto f_simplex : st.skeleton_simplex_range(3)) {
      std::clog << "vertex = (";
      for (auto vertex : st.simplex_vertex_range(f_simplex)) {
        std::clog << vertex << ",";
      }
      std::clog << ") - filtration = " << st.filtration(f_simplex);
      std::clog << " - dimension = " << st.dimension(f_simplex) << std::endl;
      if (st.dimension(f_simplex) < dimension)
        BOOST_CHECK(st.filtration(f_simplex) >= OneCriticalFiltration({1.,1}));
      else
        BOOST_CHECK(st.filtration(f_simplex) == st.inf_);
    }
  }

}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_tree_clear, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST SIMPLEX TREE CLEAR" << std::endl;
  typeST st;
  st.insert_simplex_and_subfaces({0, 1}, OneCriticalFiltration({1.5}));
  st.initialize_filtration();
  st.clear();
  BOOST_CHECK(st.num_vertices() == 0);
  BOOST_CHECK(st.num_simplices() == 0);
  BOOST_CHECK(st.upper_bound_dimension() == -1);
  BOOST_CHECK(st.dimension() == -1);
  BOOST_CHECK(boost::size(st.filtration_simplex_range()) == 0);
  typeST st_empty;
  BOOST_CHECK(st == st_empty);
  st.insert_simplex_and_subfaces({0}, OneCriticalFiltration({2.5}));
  BOOST_CHECK(boost::size(st.cofaces_simplex_range(st.find({0}), 1)) == 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multify_simplex_tree, typeST, list_of_tested_variants) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST MULTIFY FLATTEN LINEAR PROJECTION" << std::endl;
 
  typeST_STD st;
  typeST st_multi;
  st.insert_simplex_and_subfaces({1,2,3}, 1.);
  BOOST_CHECK(st.get_number_of_parameters() == 1);
  int num_parameters = 3;
  Gudhi::multi_persistence::multify(st,st_multi,num_parameters,{2.,3.}); //fills st_multi by simplices of filtration [{st.filtration(sh)}, default_values[0],default_values[1], ...]
  BOOST_CHECK(st_multi.get_number_of_parameters() == num_parameters); // num parameters is defined by multify
  BOOST_CHECK(st_multi.num_simplices() == st.num_simplices()); // simplicial complexes should be the same
  for (auto sh : st_multi.complex_simplex_range()){
    const auto& filtration = st_multi.filtration(sh);
    BOOST_CHECK(filtration == OneCriticalFiltration({1,2,3})); // Checks the filtration values
  }
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST FLATTEN" << std::endl;

  for (int dimension=0; dimension < 3; dimension++){
    st.clear();
    Gudhi::multi_persistence::flatten(st, st_multi, dimension);
    for (auto sh : st.complex_simplex_range()){
      BOOST_CHECK(st.filtration(sh) == dimension +1);
    }
  }

  // std::clog << "********************************************************************" << std::endl;
  // std::clog << "TEST LINEAR PROJECTION" << std::endl;
  // // st has already the same simplicial complex as st_multi, and the filtration values of st_multi are all {1,2,3}
  // Gudhi::multi_persistence::linear_projection(st,st_multi,{17,37,73}); // sets the filtration values of st to the dot product of st_multi.filtration and {17,37,73}.
  // for (auto sh : st.complex_simplex_range()){
  //   BOOST_CHECK(st.filtration(sh) == 1*17 + 2*37 + 3*73);
  // }
}




BOOST_AUTO_TEST_CASE(simplex_tree_multi_assign_filtration) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST MULTI FILTRATION INSERT SIMPLEX AND SUBFACES" << std::endl;
  Stree_multi st;

  typename Stree_multi::Simplex_handle sh;
  bool success = false;
  const OneCriticalFiltration multi_filt_1 = {1., 2.};
  std::tie(sh, success) = st.insert_simplex_and_subfaces({0, 1}, multi_filt_1);
  BOOST_CHECK(success);
  BOOST_CHECK(sh != st.null_simplex());
  // Only [0,1], [0] and [1] are already inserted
  const OneCriticalFiltration multi_filt_2 = {3., 2., 1.};
  std::tie(sh, success) = st.insert_simplex_and_subfaces({2, 1, 0}, multi_filt_2);
  BOOST_CHECK(success);
  BOOST_CHECK(sh != st.null_simplex());
  // Already inserted
  std::tie(sh, success) = st.insert_simplex_and_subfaces({0, 2}, {4.});
  BOOST_CHECK(!success);
  BOOST_CHECK(sh != st.null_simplex());

  // Check filtration values of an already inserted simplex
  sh = st.find({1, 2});
  BOOST_CHECK(sh != st.null_simplex());
  BOOST_CHECK(st.filtration(sh) == multi_filt_2);
  // And assign a new value
  OneCriticalFiltration const multi_filt_3 = {5.};
  st.assign_filtration(sh, multi_filt_3);

  for (auto f_simplex : st.complex_simplex_range()) {
    std::clog << "vertex = (";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << ",";
    }
    std::clog << ") - filtration = " << st.filtration(f_simplex);
    std::clog << " - dimension = " << st.dimension(f_simplex) << std::endl;
  }

  // Check all filtration values
  sh = st.find({2, 1, 0});
  BOOST_CHECK(sh != st.null_simplex());
  BOOST_CHECK(st.filtration(sh) == multi_filt_2);
  sh = st.find({1, 0});
  BOOST_CHECK(sh != st.null_simplex());
  BOOST_CHECK(st.filtration(sh) == multi_filt_1);
  sh = st.find({2, 0});
  BOOST_CHECK(sh != st.null_simplex());
  BOOST_CHECK(st.filtration(sh) == multi_filt_2);
  sh = st.find({0});
  BOOST_CHECK(sh != st.null_simplex());
  BOOST_CHECK(st.filtration(sh) == multi_filt_1);
  sh = st.find({2, 1});
  BOOST_CHECK(sh != st.null_simplex());
  BOOST_CHECK(st.filtration(sh) == multi_filt_3);
  sh = st.find({1});
  BOOST_CHECK(sh != st.null_simplex());
  BOOST_CHECK(st.filtration(sh) == multi_filt_1);
  sh = st.find({2});
  BOOST_CHECK(sh != st.null_simplex());
  BOOST_CHECK(st.filtration(sh) == multi_filt_2);

}

BOOST_AUTO_TEST_CASE(simplex_tree_multi_reset_filtration) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST RESET MULTI FILTRATION" << std::endl;
  Stree_multi st;

  st.insert_simplex_and_subfaces({2, 1, 0}, {3., 2.});
  st.insert_simplex_and_subfaces({3, 0}, {2., 3.});
  st.insert_simplex_and_subfaces({3, 4, 5}, {3., 2.});
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, {4.,4.});

  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o           */
  /*  2   0   3\X/4           */
  /*            o             */
  /*            5             */

  for (auto f_simplex : st.complex_simplex_range()) {
    std::clog << "vertex = (";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << ",";
    }
    std::clog << ") - filtration = " << st.filtration(f_simplex);
    std::clog << " - dimension = " << st.dimension(f_simplex) << std::endl;
    // Guaranteed by construction
    BOOST_CHECK(st.filtration(f_simplex) >= 2.);
  }
 
  OneCriticalFiltration new_filt = {0., 1., 2.};
  // dimension until 5 even if simplex tree is of dimension 3 to test the limits
  for(int dimension = 5; dimension >= 0; dimension --) {
    std::clog << "### reset_filtration - dimension = " << dimension << "\n";
    st.reset_filtration(new_filt, dimension);
    for (auto f_simplex : st.complex_simplex_range()) {
      std::clog << "vertex = (";
      for (auto vertex : st.simplex_vertex_range(f_simplex)) {
        std::clog << vertex << ",";
      }
      std::clog << ") - filtration = " << st.filtration(f_simplex);
      std::clog << " - dimension = " << st.dimension(f_simplex) << std::endl;
      if (st.dimension(f_simplex) >= dimension)
        BOOST_CHECK(st.filtration(f_simplex) == new_filt);
    }
  }
}

BOOST_AUTO_TEST_CASE(simplex_tree_multi_filtration_multify) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST MULTI FILTRATION MULTIFY" << std::endl;

  Simplex_tree<> st;
  st.insert_simplex_and_subfaces({2, 1, 0}, 3.);
  st.insert_simplex_and_subfaces({3, 0}, 2.);
  st.insert_simplex_and_subfaces({3, 4, 5}, 3.);
  st.insert_simplex_and_subfaces({0, 1, 6, 7}, 4.);

  /* Inserted simplex:        */
  /*    1   6                 */
  /*    o---o                 */
  /*   /X\7/                  */
  /*  o---o---o---o           */
  /*  2   0   3\X/4           */
  /*            o             */
  /*            5             */

  Stree_multi st_multi;
  OneCriticalFiltration default_multi (3, std::numeric_limits<Stree_multi::Options::value_type>::quiet_NaN());
  Gudhi::multi_persistence::multify(st, st_multi, 3, default_multi);

  for (auto f_simplex : st_multi.complex_simplex_range()) {
    std::clog << "vertex = (";
    for (auto vertex : st_multi.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << ",";
    }
    auto multi_filtration = st_multi.filtration(f_simplex);
    std::clog << ") - filtration = " << multi_filtration << std::endl;
    BOOST_CHECK(!std::isnan(multi_filtration[0]));
    BOOST_CHECK(std::isnan(multi_filtration[1]));
    BOOST_CHECK(std::isnan(multi_filtration[2]));
  }

  {
    Simplex_tree<> copy;
    Gudhi::multi_persistence::flatten(copy, st_multi, 0);
    BOOST_CHECK(st == copy);
  }

  {
    Simplex_tree<> copy;
    Gudhi::multi_persistence::flatten(copy, st_multi, 1);
    for (auto f_simplex : copy.complex_simplex_range()) {
      std::clog << "vertex = (";
      for (auto vertex : copy.simplex_vertex_range(f_simplex)) {
        std::clog << vertex << ",";
      }
      std::clog << ") - filtration = " << copy.filtration(f_simplex) << std::endl;
      BOOST_CHECK(std::isnan(copy.filtration(f_simplex)));
    }
  }

}

BOOST_AUTO_TEST_CASE(simplex_tree_multi_filtration_numeric_limits) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST MULTI FILTRATION NUMERIC LIMITS" << std::endl;

  // NaN
  auto nan_multi = std::numeric_limits<OneCriticalFiltration>::quiet_NaN();
  BOOST_CHECK(nan_multi.size() == 1);
  BOOST_CHECK(std::isnan(nan_multi[0]));
  std::clog << nan_multi << std::endl;

  // Inf
  auto inf_multi = std::numeric_limits<OneCriticalFiltration>::infinity();
  BOOST_CHECK(inf_multi.size() == 1);
  BOOST_CHECK(std::isinf(inf_multi[0]));
  std::clog << inf_multi << std::endl;
}

BOOST_AUTO_TEST_CASE(make_filtration_non_decreasing_on_multi) {
  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST MULTI FILTRATION NON DECREASING" << std::endl;
  Stree_multi st;

  st.insert_simplex_and_subfaces({2, 1, 0}, {3., 4.});
  st.insert_simplex_and_subfaces({3, 0},    {2., 3.});
  st.insert_simplex_and_subfaces({3, 4, 5}, {3., 4.});
  
  /* Inserted simplex:     */
  /*    1                  */
  /*    o                  */
  /*   /X\                 */
  /*  o---o---o---o        */
  /*  2   0   3\X/4        */
  /*            o          */
  /*            5          */


  for (auto f_simplex : st.complex_simplex_range()) {
    std::clog << "vertex = (";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << ",";
    }
    std::clog << ") - filtration = " << st.filtration(f_simplex) << std::endl;
  }
  auto filt = st.filtration(st.find({3, 0}));
  std::clog << "filtration([3,0]) = " << filt << std::endl;
  BOOST_CHECK(filt == OneCriticalFiltration({2., 3.}));

  st.set_number_of_parameters(2);
  st.make_filtration_non_decreasing();

  filt = st.filtration(st.find({3, 0}));
  std::clog << "filtration([3,0]) = " << filt << std::endl;
  BOOST_CHECK(filt == OneCriticalFiltration({3., 4.}));

  st.assign_filtration(st.find({3, 0}), OneCriticalFiltration({2.9, 4.}));
  filt = st.filtration(st.find({3, 0}));
  std::clog << "filtration([3,0]) = " << filt << std::endl;
  BOOST_CHECK(filt == OneCriticalFiltration({2.9, 4.}));

  st.make_filtration_non_decreasing();

  filt = st.filtration(st.find({3, 0}));
  std::clog << "filtration([3,0]) = " << filt << std::endl;
  BOOST_CHECK(filt == OneCriticalFiltration({3., 4.}));

  st.assign_filtration(st.find({3, 0}), OneCriticalFiltration({5., 3.99}));
  filt = st.filtration(st.find({3, 0}));
  std::clog << "filtration([3,0]) = " << filt << std::endl;
  BOOST_CHECK(filt == OneCriticalFiltration({5., 3.99}));
  
  st.make_filtration_non_decreasing();

  filt = st.filtration(st.find({3, 0}));
  std::clog << "filtration([3,0]) = " << filt << std::endl;
  BOOST_CHECK(filt == OneCriticalFiltration({5., 4.}));
}

BOOST_AUTO_TEST_CASE(make_filtration_non_decreasing_on_multi_nan_values) {
  Stree_multi st;
  
  BOOST_CHECK(std::numeric_limits<OneCriticalFiltration>::quiet_NaN().is_nan());
  BOOST_CHECK(std::numeric_limits<OneCriticalFiltration>::infinity().is_inf());

  st.insert_simplex_and_subfaces({2, 1, 0}, {1.,2.,3.});
  st.insert_simplex_and_subfaces({3, 0},    {1.,2.,3.});
  st.insert_simplex_and_subfaces({3, 4, 5}, {1.,2.,3.});

  st.assign_filtration(st.find({0}), std::numeric_limits<OneCriticalFiltration>::quiet_NaN());
  st.assign_filtration(st.find({3}), std::numeric_limits<OneCriticalFiltration>::infinity());
  
  /* Inserted simplex:     */
  /*    1                  */
  /*    o                  */
  /*   /X\                 */
  /*  o---o---o---o        */
  /*  2   0   3\X/4        */
  /*            o          */
  /*            5          */

  // Default number of parameter is 2
  BOOST_CHECK(st.get_number_of_parameters() == 2);
  
  st.set_number_of_parameters(3);
  BOOST_CHECK(st.get_number_of_parameters() == 3);

  std::clog << "SPECIFIC CASE:" << std::endl;
  std::clog << "Insertion with NaN values does not ensure the filtration values are non decreasing" << std::endl;
  st.make_filtration_non_decreasing();

  std::clog << "Check that NaN filtrations are ignored, and inf filtrations are propagated." << std::endl;
  for (auto f_simplex : st.complex_simplex_range()) {
    auto filt = st.filtration(f_simplex);
    bool contains3 = false;
    bool contains0 = false;
    int dim = -1;
    std::clog << "Simplex ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)){
      std::clog << vertex << " ";
      if (vertex == 3) contains3 = true;
      else if (vertex == 0) contains0 = true;
      dim++;
    }
    bool is_zero = dim == 0 && contains0;
    std::clog << "Filtration: " << filt << std::endl;
    if (is_zero) BOOST_CHECK(filt.is_nan());
    else if (contains3) BOOST_CHECK(filt.is_inf());
    else BOOST_CHECK(filt == OneCriticalFiltration({1.,2.,3.}));
  }
}

