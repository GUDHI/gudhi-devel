/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "vineyard"
#include <boost/test/unit_test.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/vineyard_helper.h>

#include "vy_test_utilities.h"

using ST = Gudhi::Simplex_tree<>;
using CC = Gudhi::cubical_complex::Bitmap_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> >;

BOOST_AUTO_TEST_CASE(vy_helper_simplex_tree)
{
  ST st;
  st.insert_simplex_and_subfaces({0, 1, 2}, 0.5);
  st.insert_simplex_and_subfaces({0, 2, 3}, 0.9);
  st.insert_simplex_and_subfaces({0, 2}, 0.2);

  BC bc;
  DC dc;
  FC<ST::Filtration_value> fc;

  Gudhi::vineyard::build_boundary_matrix_from_complex(st, bc, dc, fc);

  BC real_bc = {
    {},
    {},
    {},
    {},
    {5, 7, 9},
    {0, 1},
    {7, 8, 10},
    {0, 2},
    {0, 3},
    {1, 2},
    {2, 3}
  };

  DC real_dc = {0, 0, 0, 0, 2, 1, 2, 1, 1, 1, 1};

  FC<ST::Filtration_value> real_fc = {0.2, 0.5, 0.2, 0.9, 0.5, 0.5, 0.9, 0.2, 0.9, 0.5, 0.9};

  BOOST_CHECK(real_bc == bc);
  BOOST_CHECK(real_dc == dc);
  BOOST_CHECK(real_fc == fc);

  fc.clear();
  Gudhi::vineyard::build_boundary_matrix_from_complex(st, fc);
  BOOST_CHECK(real_fc == fc);
}

BOOST_AUTO_TEST_CASE(vy_helper_cubical)
{
  CC cc({2, 1}, std::vector<double>{0.5, 0.9}, true);

  BC bc;
  DC dc;
  FC<ST::Filtration_value> fc;

  Gudhi::vineyard::build_boundary_matrix_from_complex(cc, bc, dc, fc);

  BC real_bc = {
    {},
    {},
    {},
    {},
    {},
    {},
    {0, 1},
    {1, 2},
    {0, 3},
    {6, 8, 10, 13},
    {1, 4},
    {7, 10, 12, 14},
    {2, 5},
    {3, 4},
    {4, 5}
  };

  DC real_dc = {0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 1, 2, 1, 1, 1};

  FC<ST::Filtration_value> real_fc = {0.5, 0.5, 0.9, 0.5, 0.5, 0.9, 0.5, 0.9, 0.5, 0.5, 0.5, 0.9, 0.9, 0.5, 0.9};

  BOOST_CHECK(real_bc == bc);
  BOOST_CHECK(real_dc == dc);
  BOOST_CHECK(real_fc == fc);

  fc.clear();
  Gudhi::vineyard::build_boundary_matrix_from_complex(cc, fc);
  BOOST_CHECK(real_fc == fc);
}
