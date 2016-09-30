/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simple_witness_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>

#include <iostream>
#include <ctime>
#include <vector>

typedef Gudhi::Simplex_tree<> Simplex_tree;
typedef std::vector< int > typeVectorVertex;
typedef Gudhi::witness_complex::Witness_complex<Simplex_tree> WitnessComplex;

BOOST_AUTO_TEST_CASE(simple_witness_complex) {
  Simplex_tree complex;
  std::vector< typeVectorVertex > knn;

  knn.push_back({1, 0, 5, 2, 6, 3, 4});
  knn.push_back({2, 6, 4, 5, 0, 1, 3});
  knn.push_back({3, 4, 2, 1, 5, 6, 0});
  knn.push_back({4, 2, 1, 3, 5, 6, 0});
  knn.push_back({5, 1, 6, 0, 2, 3, 4});
  knn.push_back({6, 0, 5, 2, 1, 3, 4});
  knn.push_back({0, 5, 6, 1, 2, 3, 4});
  knn.push_back({2, 6, 4, 5, 3, 1, 0});
  knn.push_back({1, 2, 5, 4, 3, 6, 0});
  knn.push_back({3, 4, 0, 6, 5, 1, 2});
  knn.push_back({5, 0, 1, 3, 6, 2, 4});
  knn.push_back({5, 6, 1, 0, 2, 3, 4});
  knn.push_back({1, 6, 0, 5, 2, 3, 4});
  WitnessComplex witnessComplex(knn, 7, 7, complex);

  BOOST_CHECK(witnessComplex.is_witness_complex(knn, false));
}
