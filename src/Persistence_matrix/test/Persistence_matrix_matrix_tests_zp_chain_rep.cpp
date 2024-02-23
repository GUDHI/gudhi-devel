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
#ifdef PM_TEST_MAX_DIM
using opts = boost::mp11::mp_list<false_value_list, true_value_list>;
#else
using opts = boost::mp11::mp_list<false_value_list, false_value_list>;
#endif
#else
#ifdef PM_TEST_MAX_DIM
using opts = boost::mp11::mp_list<true_value_list, true_value_list>;
#else
using opts = boost::mp11::mp_list<true_value_list, false_value_list>;
#endif
#endif

using full_matrices = matrices_list<boost::mp11::mp_apply<opt_chain_rep_zp, opts> >;
using row_access_matrices = matrices_list<boost::mp11::mp_apply<opt_chain_rep_zp_ra, opts> >;
using removable_rows_matrices = matrices_list<boost::mp11::mp_apply<opt_chain_rep_zp_ra_r, opts> >;
using removable_columns_matrices = matrices_list<boost::mp11::mp_apply<opt_chain_rep_zp_r, opts> >;
using barcode_matrices = matrices_list<boost::mp11::mp_apply<opt_chain_rep_zp_barcode, opts> >;

// #include <boost/mpl/list.hpp>
// using opt = Chain_rep_options<Column_types::INTRUSIVE_LIST,true, false, true, true, true>;
// using mat = Matrix<opt>;
// typedef boost::mpl::list<mat> test;
// typedef boost::mpl::list<Matrix<Chain_rep_options<Column_types::INTRUSIVE_LIST,true, true, true, true, true> > > test2;

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_constructors, Matrix, full_matrices) {
	test_constructors<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_insertion, Matrix, full_matrices) {
	Matrix m1;
	m1.set_characteristic(5);

	auto orderedBoundaries = build_simple_boundary_matrix<typename Matrix::Column_type>();
	orderedBoundaries.pop_back();
	orderedBoundaries.pop_back();

	Matrix m2(orderedBoundaries, 5);

	test_chain_boundary_insertion<Matrix>(m1, m2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_access, Matrix, full_matrices) {
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);
	test_chain_access<Matrix>(m);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_row_access, Matrix, row_access_matrices) {
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);
	test_non_base_row_access<Matrix>(m);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_row_removal, Matrix, removable_rows_matrices) {
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);
	test_chain_row_removal<Matrix>(m);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_column_removal, Matrix, removable_columns_matrices) {
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);
	test_chain_maximal_simplex_removal<Matrix>(m);
}

#ifdef PM_TEST_MAX_DIM
BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_max_dimension, Matrix, full_matrices) {
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);
	test_maximal_dimension<Matrix>(m);
}
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_operation, Matrix, full_matrices) {
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);
	test_chain_operation<Matrix>(m);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_barcode, Matrix, barcode_matrices) {
	test_barcode<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_representative_cycles, Matrix, full_matrices) {
	auto columns = build_longer_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);
	test_representative_cycles<Matrix>(m);
}




