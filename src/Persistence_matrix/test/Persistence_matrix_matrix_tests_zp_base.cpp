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
// #define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
// #define BOOST_MPL_LIMIT_LIST_SIZE 10
#include <boost/test/unit_test.hpp>
// #include <boost/mpl/list.hpp>
// #include <boost/mp11.hpp>

#include "pm_matrix_tests.h"
#include "pm_matrix_tests_boost_type_lists.h"
// #include "matrix_tests_options.h"

using full_matrices = matrices_list<opt_base_zp>;
using row_access_matrices = matrices_list<opt_base_zp_ra>;
using removable_rows_matrices = matrices_list<opt_base_zp_ra_r>;
using removable_columns_matrices = matrices_list<opt_base_zp_r>;
using swap_matrices = matrices_list<opt_base_zp_swap>;

// using opt = Base_options<false,Column_types::INTRUSIVE_LIST,true,true>;
// using mat = Matrix<opt>;
// typedef boost::mpl::list<mat> test;

// std::set<std::string> test;

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_zp_constructors, Matrix, full_matrices) {
	// test.insert(typeid(Matrix()).name());
	// std::cout << test.size() << "\n";
	test_constructors<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_zp_insertion, Matrix, full_matrices) {
	test_general_insertion<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_zp_access, Matrix, full_matrices) {
	test_base_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_zp_zeroing, Matrix, full_matrices) {
	test_zeroing<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_zp_row_access, Matrix, row_access_matrices) {
	test_base_z5_row_access<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_zp_row_removal, Matrix, removable_rows_matrices) {
	test_row_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_zp_column_removal, Matrix, removable_columns_matrices) {
	test_column_removal<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_zp_operation, Matrix, full_matrices) {
	test_base_operation<Matrix>();
	test_base_cell_range_operation<Matrix>();
	test_const_operation<Matrix>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_zp_swaps, Matrix, swap_matrices) {
	test_base_swaps<Matrix>();
}

// BOOST_AUTO_TEST_CASE(Base_matrix_insertion) {
// 	boost::mp11::mp_for_each<removable_rows_matrices>([](auto M){
//         test_constructors<decltype(M)>();
//     });
// }

// template<typename T>
// void test(){
// 	boost::mp11::mp_for_each<T>( [](auto I ){
//         test_constructors<decltype(I)>();
//     });
// }

// int main(int argc, char *argv[]) {
//   test<removable_rows_matrices>();

//   return 0;
// }

