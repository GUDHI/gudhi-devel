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

using full_matrices = matrices_list<z2_ru_vine_option_list>;
using row_access_matrices = matrices_list<z2_ra_ru_vine_option_list>;
using removable_rows_matrices = matrices_list<z2_ra_r_ru_vine_option_list>;
using removable_columns_matrices = matrices_list<z2_r_ru_vine_option_list>;
using max_dim_matrices = matrices_list<z2_dim_ru_vine_option_list>;
using barcode_matrices = matrices_list<z2_barcode_ru_vine_option_list>;
using rep_matrices = matrices_list<z2_rep_ru_vine_option_list>;

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_constructors, Matrix, full_matrices) {
	test_constructors<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_insertion, Matrix, full_matrices) {
	test_boundary_insertion<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_access, Matrix, full_matrices) {
	test_boundary_access<Matrix>();
	if constexpr (Matrix::Option_list::is_indexed_by_position) test_ru_u_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_zeroing, Matrix, full_matrices) {
	test_zeroing<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_row_access, Matrix, row_access_matrices) {
	test_non_base_row_access<Matrix>();
	if constexpr (Matrix::Option_list::is_indexed_by_position) test_ru_u_row_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_row_removal, Matrix, removable_rows_matrices) {
	test_row_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_column_removal, Matrix, removable_columns_matrices) {
	test_ru_maximal_simplex_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_max_dimension, Matrix, max_dim_matrices) {
	test_maximal_dimension<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_operation, Matrix, full_matrices) {
	test_ru_operation<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_barcode, Matrix, barcode_matrices) {
	test_barcode<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine, Matrix, full_matrices) {
	auto columns = build_longer_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns);

	if constexpr (Matrix::Option_list::is_indexed_by_position) {
		test_vine_swap_with_position_index(m);
	} else {
		test_vine_swap_with_id_index(m);
	}
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_vine_representative_cycles, Matrix, rep_matrices) {
	test_representative_cycles<Matrix>();
}




