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

using full_matrices = matrices_list<opt_col_comp_zp>;
using row_access_matrices = matrices_list<opt_col_comp_zp_ra>;
using removable_rows_matrices = matrices_list<opt_col_comp_zp_ra_r>;

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_compression_matrix_zp_constructors, Matrix, full_matrices) {
	test_constructors<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_compression_matrix_zp_insertion, Matrix, full_matrices) {
	test_general_insertion<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_compression_matrix_zp_access, Matrix, full_matrices) {
	test_base_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_compression_matrix_zp_row_access, Matrix, row_access_matrices) {
	test_base_z5_row_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_compression_matrix_zp_row_removal, Matrix, removable_rows_matrices) {
	test_row_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_compression_matrix_zp_operation, Matrix, full_matrices) {
	test_base_col_comp_operation<Matrix>();
	test_base_col_comp_cell_range_operation<Matrix>();
	test_base_col_comp_const_operation<Matrix>();
}


