/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "persistence_matrix"
#include <boost/test/unit_test.hpp>

#include "pm_matrix_tests.h"
#include "pm_matrix_tests_boost_type_lists.h"

using full_matrices = matrices_list<opt_base_z2>;
using row_access_matrices = matrices_list<opt_base_z2_ra>;
using removable_rows_matrices = matrices_list<opt_base_z2_ra_r>;
using removable_columns_matrices = matrices_list<opt_base_z2_r>;
using swap_matrices = matrices_list<opt_base_z2_swap>;

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_z2_constructors, Matrix, full_matrices) { test_constructors<Matrix>(); }

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_z2_insertion, Matrix, full_matrices) { test_general_insertion<Matrix>(); }

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_z2_access, Matrix, full_matrices) { test_base_access<Matrix>(); }

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_z2_zeroing, Matrix, full_matrices) { test_zeroing<Matrix>(); }

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_z2_row_access, Matrix, row_access_matrices) {
  test_base_z2_row_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_z2_row_removal, Matrix, removable_rows_matrices) {
  test_row_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_z2_column_removal, Matrix, removable_columns_matrices) {
  test_column_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_z2_operation, Matrix, full_matrices) {
  test_base_operation<Matrix>();
  test_base_entry_range_operation<Matrix>();
  test_const_operation<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_z2_swaps, Matrix, swap_matrices) { test_base_swaps<Matrix>(); }
