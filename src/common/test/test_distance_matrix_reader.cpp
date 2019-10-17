/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/reader_utils.h>

#include <iostream>
#include <string>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "distance_matrix_reader"
#include <boost/test/unit_test.hpp>

using Distance_matrix = std::vector<std::vector<double>>;

BOOST_AUTO_TEST_CASE( lower_triangular_distance_matrix )
{
  Distance_matrix from_lower_triangular;
  // Read lower_triangular_distance_matrix.csv file where the separator is a ','
  from_lower_triangular = Gudhi::read_lower_triangular_matrix_from_csv_file<double>("lower_triangular_distance_matrix.csv",
                                                                             ',');
  for (auto& i : from_lower_triangular) {
    for (auto j : i) {
      std::cout << j << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "from_lower_triangular size = " << from_lower_triangular.size() << std::endl;
  BOOST_CHECK(from_lower_triangular.size() == 5);

  for (std::size_t i = 0; i < from_lower_triangular.size(); i++) {
    std::cout << "from_lower_triangular[" << i << "] size = " << from_lower_triangular[i].size() << std::endl;
    BOOST_CHECK(from_lower_triangular[i].size() == i);
  }
  std::vector<double> expected = {1};
  BOOST_CHECK(from_lower_triangular[1] == expected);
  
  expected = {2,3};
  BOOST_CHECK(from_lower_triangular[2] == expected);
  
  expected = {4,5,6};
  BOOST_CHECK(from_lower_triangular[3] == expected);
  
  expected = {7,8,9,10};
  BOOST_CHECK(from_lower_triangular[4] == expected);
  
}

BOOST_AUTO_TEST_CASE( full_square_distance_matrix )
{
  Distance_matrix from_full_square;
  // Read full_square_distance_matrix.csv file where the separator is the default one ';'
  from_full_square = Gudhi::read_lower_triangular_matrix_from_csv_file<double>("full_square_distance_matrix.csv");
  for (auto& i : from_full_square) {
    for (auto j : i) {
      std::cout << j << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "from_full_square size = " << from_full_square.size() << std::endl;
  BOOST_CHECK(from_full_square.size() == 7);
  for (std::size_t i = 0; i < from_full_square.size(); i++) {
    std::cout << "from_full_square[" << i << "] size = " << from_full_square[i].size() << std::endl;
    BOOST_CHECK(from_full_square[i].size() == i);
  }  
}
