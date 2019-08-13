/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2017 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/reader_utils.h>

#include <iostream>
#include <vector>
#include <utility>  // for pair
#include <tuple>
#include <limits>  // for inf
#include <map>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "persistence_intervals_reader"
#include <boost/test/unit_test.hpp>

using Persistence_intervals_by_dimension = std::map<int, std::vector<std::pair<double, double>>>;
using Persistence_intervals = std::vector<std::pair<double, double>>;
// Test files with only 2 parameters (persistence birth and death) per line in file
BOOST_AUTO_TEST_CASE( persistence_intervals_without_dimension )
{
  Persistence_intervals_by_dimension expected_intervals_by_dimension;
  expected_intervals_by_dimension[-1].push_back(std::make_pair(2.7, 3.7));
  expected_intervals_by_dimension[-1].push_back(std::make_pair(9.6, 14.));
  expected_intervals_by_dimension[-1].push_back(std::make_pair(34.2, 34.974));
  expected_intervals_by_dimension[-1].push_back(std::make_pair(3., std::numeric_limits<double>::infinity()));

  Persistence_intervals_by_dimension persistence_intervals_by_dimension =
      Gudhi::read_persistence_intervals_grouped_by_dimension("persistence_intervals_without_dimension.pers");

  std::cout << "\nread_persistence_intervals_grouped_by_dimension - expected\n";
  for (auto map_iter : expected_intervals_by_dimension) {
    std::cout << "key=" << map_iter.first;
    for (auto vec_iter : map_iter.second)
      std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";
  }

  std::cout << "\nread_persistence_intervals_grouped_by_dimension - read\n";
  for (auto map_iter : persistence_intervals_by_dimension) {
    std::cout << "key=" << map_iter.first;
    for (auto vec_iter : map_iter.second)
      std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";
  }

  BOOST_CHECK(persistence_intervals_by_dimension == expected_intervals_by_dimension);

  Persistence_intervals expected_intervals_in_dimension;
  expected_intervals_in_dimension.push_back(std::make_pair(2.7, 3.7));
  expected_intervals_in_dimension.push_back(std::make_pair(9.6, 14.));
  expected_intervals_in_dimension.push_back(std::make_pair(34.2, 34.974));
  expected_intervals_in_dimension.push_back(std::make_pair(3., std::numeric_limits<double>::infinity()));

  Persistence_intervals persistence_intervals_in_dimension =
      Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_without_dimension.pers");

  std::cout << "\nread_persistence_intervals_in_dimension - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension - read\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_without_dimension.pers", 0);
  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_without_dimension.pers", 1);
  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_without_dimension.pers", 2);
  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_without_dimension.pers", 3);
  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

}
// Test files with 3 parameters (dimension birth death) per line in file
BOOST_AUTO_TEST_CASE( persistence_intervals_with_dimension )
{
  Persistence_intervals_by_dimension expected_intervals_by_dimension;
  expected_intervals_by_dimension[0].push_back(std::make_pair(2.7, 3.7));
  expected_intervals_by_dimension[1].push_back(std::make_pair(9.6, 14.));
  expected_intervals_by_dimension[3].push_back(std::make_pair(34.2, 34.974));
  expected_intervals_by_dimension[1].push_back(std::make_pair(3., std::numeric_limits<double>::infinity()));

  Persistence_intervals_by_dimension persistence_intervals_by_dimension =
      Gudhi::read_persistence_intervals_grouped_by_dimension("persistence_intervals_with_dimension.pers");

  std::cout << "\nread_persistence_intervals_grouped_by_dimension - expected\n";
  for (auto map_iter : expected_intervals_by_dimension) {
    std::cout << "key=" << map_iter.first;
    for (auto vec_iter : map_iter.second)
      std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";
  }

  std::cout << "\nread_persistence_intervals_grouped_by_dimension - read\n";
  for (auto map_iter : persistence_intervals_by_dimension) {
    std::cout << "key=" << map_iter.first;
    for (auto vec_iter : map_iter.second)
      std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";
  }

  BOOST_CHECK(persistence_intervals_by_dimension == expected_intervals_by_dimension);

  Persistence_intervals expected_intervals_in_dimension;
  expected_intervals_in_dimension.push_back(std::make_pair(2.7, 3.7));
  expected_intervals_in_dimension.push_back(std::make_pair(9.6, 14.));
  expected_intervals_in_dimension.push_back(std::make_pair(34.2, 34.974));
  expected_intervals_in_dimension.push_back(std::make_pair(3., std::numeric_limits<double>::infinity()));

  Persistence_intervals persistence_intervals_in_dimension =
      Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_dimension.pers");

  std::cout << "\nread_persistence_intervals_in_dimension - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  expected_intervals_in_dimension.push_back(std::make_pair(2.7, 3.7));
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_dimension.pers", 0);

  std::cout << "\nread_persistence_intervals_in_dimension 0 - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension 0 - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  expected_intervals_in_dimension.push_back(std::make_pair(9.6, 14.));
  expected_intervals_in_dimension.push_back(std::make_pair(3., std::numeric_limits<double>::infinity()));
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_dimension.pers", 1);

  std::cout << "\nread_persistence_intervals_in_dimension 1 - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension 1 - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_dimension.pers", 2);

  std::cout << "\nread_persistence_intervals_in_dimension 2 - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension 2 - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  expected_intervals_in_dimension.push_back(std::make_pair(34.2, 34.974));
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_dimension.pers", 3);

  std::cout << "\nread_persistence_intervals_in_dimension 3 - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension 3 - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

}

// Test files with 4 parameters (field dimension birth death) per line in file
BOOST_AUTO_TEST_CASE( persistence_intervals_with_field )
{
  Persistence_intervals_by_dimension expected_intervals_by_dimension;
  expected_intervals_by_dimension[0].push_back(std::make_pair(2.7, 3.7));
  expected_intervals_by_dimension[1].push_back(std::make_pair(9.6, 14.));
  expected_intervals_by_dimension[3].push_back(std::make_pair(34.2, 34.974));
  expected_intervals_by_dimension[1].push_back(std::make_pair(3., std::numeric_limits<double>::infinity()));

  Persistence_intervals_by_dimension persistence_intervals_by_dimension =
      Gudhi::read_persistence_intervals_grouped_by_dimension("persistence_intervals_with_field.pers");

  std::cout << "\nread_persistence_intervals_grouped_by_dimension - expected\n";
  for (auto map_iter : expected_intervals_by_dimension) {
    std::cout << "key=" << map_iter.first;
    for (auto vec_iter : map_iter.second)
      std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";
  }

  std::cout << "\nread_persistence_intervals_grouped_by_dimension - read\n";
  for (auto map_iter : persistence_intervals_by_dimension) {
    std::cout << "key=" << map_iter.first;
    for (auto vec_iter : map_iter.second)
      std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";
  }

  BOOST_CHECK(persistence_intervals_by_dimension == expected_intervals_by_dimension);

  Persistence_intervals expected_intervals_in_dimension;
  expected_intervals_in_dimension.push_back(std::make_pair(2.7, 3.7));
  expected_intervals_in_dimension.push_back(std::make_pair(9.6, 14.));
  expected_intervals_in_dimension.push_back(std::make_pair(34.2, 34.974));
  expected_intervals_in_dimension.push_back(std::make_pair(3., std::numeric_limits<double>::infinity()));

  Persistence_intervals persistence_intervals_in_dimension =
      Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_field.pers");

  std::cout << "\nread_persistence_intervals_in_dimension - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  expected_intervals_in_dimension.push_back(std::make_pair(2.7, 3.7));
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_field.pers", 0);

  std::cout << "\nread_persistence_intervals_in_dimension 0 - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension 0 - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  expected_intervals_in_dimension.push_back(std::make_pair(9.6, 14.));
  expected_intervals_in_dimension.push_back(std::make_pair(3., std::numeric_limits<double>::infinity()));
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_field.pers", 1);

  std::cout << "\nread_persistence_intervals_in_dimension 1 - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension 1 - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_field.pers", 2);

  std::cout << "\nread_persistence_intervals_in_dimension 2 - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension 2 - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

  expected_intervals_in_dimension.clear();
  expected_intervals_in_dimension.push_back(std::make_pair(34.2, 34.974));
  persistence_intervals_in_dimension =
        Gudhi::read_persistence_intervals_in_dimension("persistence_intervals_with_field.pers", 3);

  std::cout << "\nread_persistence_intervals_in_dimension 3 - expected\n";
  for (auto vec_iter : expected_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  std::cout << "\nread_persistence_intervals_in_dimension 3 - read\n";
  for (auto vec_iter : persistence_intervals_in_dimension)
    std::cout << " [" << vec_iter.first << " ," << vec_iter.second << "] ";

  BOOST_CHECK(persistence_intervals_in_dimension == expected_intervals_in_dimension);

}
