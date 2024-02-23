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
#ifdef PM_TEST_BARCODE
#ifdef PM_TEST_MAX_DIM
using opts = boost::mp11::mp_list<false_value_list, true_value_list, true_value_list>;
#else
using opts = boost::mp11::mp_list<false_value_list, true_value_list, false_value_list>;
#endif
#else
#ifdef PM_TEST_MAX_DIM
using opts = boost::mp11::mp_list<false_value_list, false_value_list, true_value_list>;
#else
using opts = boost::mp11::mp_list<false_value_list, false_value_list, false_value_list>;
#endif
#endif
#else
#ifdef PM_TEST_BARCODE
#ifdef PM_TEST_MAX_DIM
using opts = boost::mp11::mp_list<true_value_list, true_value_list, true_value_list>;
#else
using opts = boost::mp11::mp_list<true_value_list, true_value_list, false_value_list>;
#endif
#else
#ifdef PM_TEST_MAX_DIM
using opts = boost::mp11::mp_list<true_value_list, false_value_list, true_value_list>;
#else
using opts = boost::mp11::mp_list<true_value_list, false_value_list, false_value_list>;
#endif
#endif
#endif

using full_matrices = matrices_list<boost::mp11::mp_apply<opt_ru_rep_z2, opts> >;
using row_access_matrices = matrices_list<boost::mp11::mp_apply<opt_ru_rep_z2_ra, opts> >;
using removable_rows_matrices = matrices_list<boost::mp11::mp_apply<opt_ru_rep_z2_ra_r, opts> >;
using removable_columns_matrices = matrices_list<boost::mp11::mp_apply<opt_ru_rep_z2_r, opts> >;

// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_constructors, Matrix, full_matrices) {
// 	test_constructors<Matrix>();
// }

// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_insertion, Matrix, full_matrices) {
// 	test_boundary_insertion<Matrix>();
// }

// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_zeroing, Matrix, full_matrices) {
// 	test_zeroing<Matrix>();
// }

// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_row_removal, Matrix, removable_rows_matrices) {
// 	test_row_removal<Matrix>();
// }

// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_column_removal, Matrix, removable_columns_matrices) {
// 	test_ru_maximal_simplex_removal<Matrix>();
// }

// #ifdef PM_TEST_MAX_DIM
// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_max_dimension, Matrix, full_matrices) {
// 	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
// 	Matrix m(columns);
// 	test_maximal_dimension<Matrix>(m);
// }
// #endif

BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_operation, Matrix, full_matrices) {
	test_ru_operation<Matrix>();
}

// #ifdef PM_TEST_BARCODE
// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_barcode, Matrix, full_matrices) {
// 	test_barcode<Matrix>();
// }
// #endif

// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_representative_cycles, Matrix, full_matrices) {
// 	auto columns = build_longer_boundary_matrix<typename Matrix::Column_type>();
// 	Matrix m(columns);
// 	test_representative_cycles<Matrix>(m);
// }

// #ifdef PM_TEST_ID_IDX
// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_access, Matrix, full_matrices) {
// 	test_boundary_access<Matrix>();
// }

// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_row_access, Matrix, row_access_matrices) {
// 	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
// 	Matrix m(columns);
// 	test_non_base_row_access<Matrix>(m);
// }
// #else
// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_access, Matrix, full_matrices) {
// 	test_boundary_access<Matrix>();
// 	test_ru_u_access<Matrix>();
// }

// BOOST_AUTO_TEST_CASE_TEMPLATE(RU_matrix_z2_rep_row_access, Matrix, row_access_matrices) {
// 	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
// 	Matrix m(columns);
// 	test_non_base_row_access<Matrix>(m);
// 	test_ru_u_row_access<Matrix>();
// }
// #endif




