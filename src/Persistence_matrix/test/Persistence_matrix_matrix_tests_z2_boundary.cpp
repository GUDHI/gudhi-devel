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

#ifdef PM_TEST_ID_IDX
using full_matrices = matrices_list<opt_boundary_z2<false_value_list> >;
using row_access_matrices = matrices_list<opt_boundary_z2_ra<false_value_list> >;
using removable_rows_matrices = matrices_list<opt_boundary_z2_ra_r<false_value_list> >;
using removable_columns_matrices = matrices_list<opt_boundary_z2_r<false_value_list> >;
using max_dim_matrices = matrices_list<opt_boundary_z2_dim<false_value_list> >;
using barcode_matrices = matrices_list<opt_boundary_z2_barcode<false_value_list> >;
using swap_matrices = matrices_list<opt_boundary_z2_swap<false_value_list> >;
#else
using full_matrices = matrices_list<opt_boundary_z2<true_value_list> >;
using row_access_matrices = matrices_list<opt_boundary_z2_ra<true_value_list> >;
using removable_rows_matrices = matrices_list<opt_boundary_z2_ra_r<true_value_list> >;
using removable_columns_matrices = matrices_list<opt_boundary_z2_r<true_value_list> >;
using max_dim_matrices = matrices_list<opt_boundary_z2_dim<true_value_list> >;
using barcode_matrices = matrices_list<opt_boundary_z2_barcode<true_value_list> >;
using swap_matrices = matrices_list<opt_boundary_z2_swap<true_value_list> >;
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_constructors, Matrix, full_matrices) {
	test_constructors<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_insertion, Matrix, full_matrices) {
	test_boundary_insertion<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_access, Matrix, full_matrices) {
	test_boundary_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_zeroing, Matrix, full_matrices) {
	test_zeroing<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_row_access, Matrix, row_access_matrices) {
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns);
	test_non_base_row_access<Matrix>(m);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_row_removal, Matrix, removable_rows_matrices) {
	test_row_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_column_removal, Matrix, removable_columns_matrices) {
	test_boundary_maximal_simplex_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_max_dimension, Matrix, max_dim_matrices) {
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns);
	test_maximal_dimension<Matrix>(m);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_operation, Matrix, full_matrices) {
	test_base_operation<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_barcode, Matrix, barcode_matrices) {
	test_barcode<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_z2_swaps, Matrix, swap_matrices) {
	test_base_swaps<Matrix>();
}



