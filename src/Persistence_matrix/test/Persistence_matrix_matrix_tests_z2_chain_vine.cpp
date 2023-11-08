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

using row_access_matrices = matrices_list<z2_ra_chain_vine_option_list>;
using removable_rows_matrices = matrices_list<z2_ra_r_chain_vine_option_list>;
using removable_columns_matrices = matrices_list<z2_r_chain_vine_option_list>;
using max_dim_matrices = matrices_list<z2_dim_chain_vine_option_list>;
using barcode_matrices = matrices_list<z2_barcode_chain_vine_option_list>;
using no_barcode_matrices = matrices_list<z2_no_barcode_chain_vine_option_list>;
using rep_matrices = matrices_list<z2_rep_chain_vine_option_list>;

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_constructors, Matrix, barcode_matrices) {
	test_constructors<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_no_barcode_constructors, Matrix, no_barcode_matrices) {
	test_chain_constructors<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_insertion, Matrix, barcode_matrices) {
	test_chain_boundary_insertion<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_access, Matrix, barcode_matrices) {
	test_chain_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_row_access, Matrix, row_access_matrices) {
	test_non_base_row_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_row_removal, Matrix, removable_rows_matrices) {
	test_chain_row_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_column_removal, Matrix, removable_columns_matrices) {
	test_chain_maximal_simplex_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_max_dimension, Matrix, max_dim_matrices) {
	test_maximal_dimension<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_operation, Matrix, barcode_matrices) {
	test_chain_operation<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_barcode, Matrix, barcode_matrices) {
	test_barcode<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine, Matrix, barcode_matrices) {
	auto columns = build_longer_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns);

	if constexpr (Matrix::Option_list::is_indexed_by_position) {
		test_vine_swap_with_position_index(m);
	} else {
		test_vine_swap_with_id_index(m);
	}
}

//only works because of the example used
bool birth_comparator(unsigned int columnIndex1, unsigned int columnIndex2) {
	if (columnIndex1 == 0) return false;
	if (columnIndex1 == 3) return true;
	return false;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_without_barcode, Matrix, no_barcode_matrices) {
	auto columns = build_longer_boundary_matrix<typename Matrix::Column_type>();

	Matrix m(columns,birth_comparator,Gudhi::persistence_matrix::_no_G_death_comparator);

	if constexpr (Matrix::Option_list::is_indexed_by_position) {
		test_vine_swap_with_position_index(m);
	} else {
		test_vine_swap_with_id_index(m);
	}
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_z2_vine_representative_cycles, Matrix, rep_matrices) {
	test_representative_cycles<Matrix>();
}




