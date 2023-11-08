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

using full_matrices = matrices_list<zp_chain_rep_option_list>;
using row_access_matrices = matrices_list<zp_ra_chain_rep_option_list>;
using removable_rows_matrices = matrices_list<zp_ra_r_chain_rep_option_list>;
using removable_columns_matrices = matrices_list<zp_r_chain_rep_option_list>;
using max_dim_matrices = matrices_list<zp_dim_chain_rep_option_list>;
using barcode_matrices = matrices_list<zp_barcode_chain_rep_option_list>;

// #include <boost/mpl/list.hpp>
// using opt = Chain_rep_options<Column_types::INTRUSIVE_LIST,true, false, true, true, true>;
// using mat = Matrix<opt>;
// typedef boost::mpl::list<mat> test;
// typedef boost::mpl::list<Matrix<Chain_rep_options<Column_types::INTRUSIVE_LIST,true, true, true, true, true> > > test2;

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_constructors, Matrix, full_matrices) {
	test_constructors<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_insertion, Matrix, full_matrices) {
	test_chain_boundary_insertion<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_access, Matrix, full_matrices) {
	test_chain_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_row_access, Matrix, row_access_matrices) {
	test_non_base_row_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_row_removal, Matrix, removable_rows_matrices) {
	test_chain_row_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_column_removal, Matrix, removable_columns_matrices) {
	test_chain_maximal_simplex_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_max_dimension, Matrix, max_dim_matrices) {
	test_maximal_dimension<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_operation, Matrix, full_matrices) {
	test_chain_operation<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_barcode, Matrix, barcode_matrices) {
	test_barcode<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_zp_rep_representative_cycles, Matrix, full_matrices) {
	test_representative_cycles<Matrix>();
}




