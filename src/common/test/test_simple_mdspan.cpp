/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_mdspan"
#include <boost/test/unit_test.hpp>

#include <gudhi/Simple_mdspan.h>

using Gudhi::Simple_mdspan;

BOOST_AUTO_TEST_CASE(test_simple_mdspan_access)
{
  std::vector<int> v{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

  auto ms2 = Simple_mdspan<int>(v.data(), 2, 6);
  auto ms3 = Simple_mdspan(v.data(), 2, 3, 2);
  auto ms4 = Simple_mdspan(v.data(), {1, 12});

  BOOST_CHECK_EQUAL((ms2[{0, 0}]), 1);
  BOOST_CHECK_EQUAL((ms2[{0, 1}]), 2);
  BOOST_CHECK_EQUAL((ms2[{0, 2}]), 3);
  BOOST_CHECK_EQUAL((ms2[{0, 3}]), 4);
  BOOST_CHECK_EQUAL((ms2[{0, 4}]), 5);
  BOOST_CHECK_EQUAL((ms2[{0, 5}]), 6);
  BOOST_CHECK_EQUAL((ms2[{1, 0}]), 7);
  BOOST_CHECK_EQUAL((ms2[{1, 1}]), 8);
  BOOST_CHECK_EQUAL((ms2[{1, 2}]), 9);
  BOOST_CHECK_EQUAL((ms2[{1, 3}]), 10);
  BOOST_CHECK_EQUAL((ms2[{1, 4}]), 11);
  BOOST_CHECK_EQUAL((ms2[{1, 5}]), 12);

  BOOST_CHECK_EQUAL((ms3[{0, 0, 0}]), 1);
  BOOST_CHECK_EQUAL((ms3[{0, 0, 1}]), 2);
  BOOST_CHECK_EQUAL((ms3[{0, 1, 0}]), 3);
  BOOST_CHECK_EQUAL((ms3[{0, 1, 1}]), 4);
  BOOST_CHECK_EQUAL((ms3[{0, 2, 0}]), 5);
  BOOST_CHECK_EQUAL((ms3[{0, 2, 1}]), 6);
  BOOST_CHECK_EQUAL((ms3[{1, 0, 0}]), 7);
  BOOST_CHECK_EQUAL((ms3[{1, 0, 1}]), 8);
  BOOST_CHECK_EQUAL((ms3[{1, 1, 0}]), 9);
  BOOST_CHECK_EQUAL((ms3[{1, 1, 1}]), 10);
  BOOST_CHECK_EQUAL((ms3[{1, 2, 0}]), 11);
  BOOST_CHECK_EQUAL((ms3[{1, 2, 1}]), 12);

  BOOST_CHECK_EQUAL((ms4(0, 0)), 1);
  BOOST_CHECK_EQUAL((ms4(0, 1)), 2);
  BOOST_CHECK_EQUAL((ms4(0, 2)), 3);
  BOOST_CHECK_EQUAL((ms4(0, 3)), 4);
  BOOST_CHECK_EQUAL((ms4(0, 4)), 5);
  BOOST_CHECK_EQUAL((ms4(0, 5)), 6);
  BOOST_CHECK_EQUAL((ms4(0, 6)), 7);
  BOOST_CHECK_EQUAL((ms4(0, 7)), 8);
  BOOST_CHECK_EQUAL((ms4(0, 8)), 9);
  BOOST_CHECK_EQUAL((ms4(0, 9)), 10);
  BOOST_CHECK_EQUAL((ms4(0, 10)), 11);
  BOOST_CHECK_EQUAL((ms4(0, 11)), 12);

  for (std::size_t i = 0; i != ms2.extent(0); i++) {
    for (std::size_t j = 0; j != ms2.extent(1); j++) {
      ms2[{i, j}] = i * 1000 + j;
    }
  }

  BOOST_CHECK_EQUAL((ms2[{0, 0}]), 0);
  BOOST_CHECK_EQUAL((ms2[{0, 1}]), 1);
  BOOST_CHECK_EQUAL((ms2[{0, 2}]), 2);
  BOOST_CHECK_EQUAL((ms2[{0, 3}]), 3);
  BOOST_CHECK_EQUAL((ms2[{0, 4}]), 4);
  BOOST_CHECK_EQUAL((ms2[{0, 5}]), 5);
  BOOST_CHECK_EQUAL((ms2[{1, 0}]), 1000);
  BOOST_CHECK_EQUAL((ms2[{1, 1}]), 1001);
  BOOST_CHECK_EQUAL((ms2[{1, 2}]), 1002);
  BOOST_CHECK_EQUAL((ms2[{1, 3}]), 1003);
  BOOST_CHECK_EQUAL((ms2[{1, 4}]), 1004);
  BOOST_CHECK_EQUAL((ms2[{1, 5}]), 1005);

  BOOST_CHECK_EQUAL((ms3[{0, 0, 0}]), 0);
  BOOST_CHECK_EQUAL((ms3[{0, 0, 1}]), 1);
  BOOST_CHECK_EQUAL((ms3[{0, 1, 0}]), 2);
  BOOST_CHECK_EQUAL((ms3[{0, 1, 1}]), 3);
  BOOST_CHECK_EQUAL((ms3[{0, 2, 0}]), 4);
  BOOST_CHECK_EQUAL((ms3[{0, 2, 1}]), 5);
  BOOST_CHECK_EQUAL((ms3[{1, 0, 0}]), 1000);
  BOOST_CHECK_EQUAL((ms3[{1, 0, 1}]), 1001);
  BOOST_CHECK_EQUAL((ms3[{1, 1, 0}]), 1002);
  BOOST_CHECK_EQUAL((ms3[{1, 1, 1}]), 1003);
  BOOST_CHECK_EQUAL((ms3[{1, 2, 0}]), 1004);
  BOOST_CHECK_EQUAL((ms3[{1, 2, 1}]), 1005);

  BOOST_CHECK_EQUAL((ms4(0, 0)), 0);
  BOOST_CHECK_EQUAL((ms4(0, 1)), 1);
  BOOST_CHECK_EQUAL((ms4(0, 2)), 2);
  BOOST_CHECK_EQUAL((ms4(0, 3)), 3);
  BOOST_CHECK_EQUAL((ms4(0, 4)), 4);
  BOOST_CHECK_EQUAL((ms4(0, 5)), 5);
  BOOST_CHECK_EQUAL((ms4(0, 6)), 1000);
  BOOST_CHECK_EQUAL((ms4(0, 7)), 1001);
  BOOST_CHECK_EQUAL((ms4(0, 8)), 1002);
  BOOST_CHECK_EQUAL((ms4(0, 9)), 1003);
  BOOST_CHECK_EQUAL((ms4(0, 10)), 1004);
  BOOST_CHECK_EQUAL((ms4(0, 11)), 1005);
}

BOOST_AUTO_TEST_CASE(test_simple_mdspan_properties)
{
  std::vector<int> v{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

  auto ms2 = Simple_mdspan<int>(v.data(), 2, 6);
  auto ms3 = Simple_mdspan(v.data(), 2, 3, 2);
  auto ms4 = Simple_mdspan(v.data(), {1, 12});

  BOOST_CHECK_EQUAL(ms2.rank(), 2);
  BOOST_CHECK_EQUAL(ms3.rank(), 3);
  BOOST_CHECK_EQUAL(ms4.rank(), 2);
  BOOST_CHECK_EQUAL(ms2.rank_dynamic(), 2);
  BOOST_CHECK_EQUAL(ms3.rank_dynamic(), 3);
  BOOST_CHECK_EQUAL(ms4.rank_dynamic(), 2);

  BOOST_CHECK_EQUAL(ms2.extent(0), 2);
  BOOST_CHECK_EQUAL(ms2.extent(1), 6);
  BOOST_CHECK_EQUAL(ms3.extent(0), 2);
  BOOST_CHECK_EQUAL(ms3.extent(1), 3);
  BOOST_CHECK_EQUAL(ms3.extent(2), 2);
  BOOST_CHECK_EQUAL(ms4.extent(0), 1);
  BOOST_CHECK_EQUAL(ms4.extent(1), 12);

  BOOST_CHECK_EQUAL(ms2.size(), 12);
  BOOST_CHECK_EQUAL(ms3.size(), 12);
  BOOST_CHECK_EQUAL(ms4.size(), 12);

  BOOST_CHECK_EQUAL(ms2.empty(), false);
  BOOST_CHECK_EQUAL(ms3.empty(), false);
  BOOST_CHECK_EQUAL(ms4.empty(), false);

  BOOST_CHECK_EQUAL(ms2.stride(0), 6);
  BOOST_CHECK_EQUAL(ms2.stride(1), 1);
  BOOST_CHECK_EQUAL(ms3.stride(0), 6);
  BOOST_CHECK_EQUAL(ms3.stride(1), 2);
  BOOST_CHECK_EQUAL(ms3.stride(2), 1);
  BOOST_CHECK_EQUAL(ms4.stride(0), 12);
  BOOST_CHECK_EQUAL(ms4.stride(1), 1);
}
