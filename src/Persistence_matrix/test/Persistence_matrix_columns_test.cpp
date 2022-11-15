/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "persistence_matrix"
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_LIST_SIZE 10
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>
#include <random>
#include <vector>
#include <set>
#include <utility>

#include "gudhi/utilities/Zp_field.h"
#include "gudhi/utilities/utilities.h"
#include "gudhi/matrix.h"
#include "gudhi/options.h"
#include "gudhi/column_types/list_column.h"
#include "gudhi/column_types/set_column.h"
#include "gudhi/column_types/unordered_set_column.h"
#include "gudhi/column_types/vector_column.h"
#include "gudhi/column_types/reduced_cell_list_column_with_row.h"
#include "gudhi/column_types/reduced_cell_set_column_with_row.h"
#include "gudhi/column_types/z2_heap_column.h"
#include "gudhi/column_types/z2_list_column.h"
#include "gudhi/column_types/z2_set_column.h"
#include "gudhi/column_types/z2_unordered_set_column.h"
#include "gudhi/column_types/z2_vector_column.h"
#include "gudhi/column_types/z2_reduced_cell_list_column_with_row.h"
#include "gudhi/column_types/z2_reduced_cell_set_column_with_row.h"
#include "gudhi/column_types/column_pairing.h"

using Gudhi::persistence_matrix::Matrix;
using Gudhi::persistence_matrix::Zp_field_element;
using Gudhi::persistence_matrix::List_column;
using Gudhi::persistence_matrix::Set_column;
using Gudhi::persistence_matrix::Unordered_set_column;
using Gudhi::persistence_matrix::Vector_column;
using Gudhi::persistence_matrix::Reduced_cell_list_column_with_row;
using Gudhi::persistence_matrix::Reduced_cell_set_column_with_row;
using Gudhi::persistence_matrix::Z2_heap_column;
using Gudhi::persistence_matrix::Z2_list_column;
using Gudhi::persistence_matrix::Z2_set_column;
using Gudhi::persistence_matrix::Z2_unordered_set_column;
using Gudhi::persistence_matrix::Z2_vector_column;
using Gudhi::persistence_matrix::Z2_reduced_cell_list_column_with_row;
using Gudhi::persistence_matrix::Z2_reduced_cell_set_column_with_row;
using Gudhi::persistence_matrix::Column_pairing;
using Gudhi::persistence_matrix::Column_types;
using Gudhi::persistence_matrix::Default_options;
using Gudhi::persistence_matrix::dimension_type;

using Z5 = Zp_field_element<5>;
using Z2 = Zp_field_element<2>;

struct set_5_test_options_with_pairing : Default_options<Z5, Column_types::SET, false, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_row_access = true;
};

struct set_5_test_options_without_pairing : Default_options<Z5, Column_types::SET, false, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_row_access = true;
};

struct list_5_test_options_with_pairing : Default_options<Z5, Column_types::LIST, false, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_row_access = true;
};

struct list_5_test_options_without_pairing : Default_options<Z5, Column_types::LIST, false, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_row_access = true;
};

struct set_2_test_options_with_pairing : Default_options<Z2, Column_types::SET, false, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_row_access = true;
};

struct set_2_test_options_without_pairing : Default_options<Z2, Column_types::SET, false, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_row_access = true;
};

struct list_2_test_options_with_pairing : Default_options<Z2, Column_types::LIST, false, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_row_access = true;
};

struct list_2_test_options_without_pairing : Default_options<Z2, Column_types::LIST, false, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_row_access = true;
};

using Set_5_matrix_with_pairing = Matrix<set_5_test_options_with_pairing>;
using Set_5_matrix_without_pairing = Matrix<set_5_test_options_without_pairing>;
using List_5_matrix_with_pairing = Matrix<list_5_test_options_with_pairing>;
using List_5_matrix_without_pairing = Matrix<list_5_test_options_without_pairing>;
using Dummy_pairing_5 = Set_5_matrix_without_pairing::Dummy_column_pairing;
using Set_2_matrix_with_pairing = Matrix<set_2_test_options_with_pairing>;
using Set_2_matrix_without_pairing = Matrix<set_2_test_options_without_pairing>;
using List_2_matrix_with_pairing = Matrix<list_2_test_options_with_pairing>;
using List_2_matrix_without_pairing = Matrix<list_2_test_options_without_pairing>;
using Dummy_pairing_2 = Set_2_matrix_without_pairing::Dummy_column_pairing;

typedef boost::mpl::list<Reduced_cell_list_column_with_row<List_5_matrix_with_pairing, Z5, Column_pairing>,
							Reduced_cell_list_column_with_row<List_5_matrix_without_pairing, Z5, List_5_matrix_without_pairing::Dummy_column_pairing>,
							Reduced_cell_set_column_with_row<Set_5_matrix_with_pairing, Z5, Column_pairing>,
							Reduced_cell_set_column_with_row<Set_5_matrix_without_pairing, Z5, Dummy_pairing_5>
						> list_of_5_columns_with_row;

typedef boost::mpl::list<Z2_reduced_cell_list_column_with_row<List_2_matrix_with_pairing, Column_pairing>,
							Z2_reduced_cell_list_column_with_row<List_2_matrix_without_pairing, List_2_matrix_without_pairing::Dummy_column_pairing>,
							Z2_reduced_cell_set_column_with_row<Set_2_matrix_with_pairing, Column_pairing>,
							Z2_reduced_cell_set_column_with_row<Set_2_matrix_without_pairing, Dummy_pairing_2>
						> list_of_2_columns_with_row;

typedef boost::mpl::list<List_column<Z5, Column_pairing>,
							List_column<Z5, Dummy_pairing_5>,
							Set_column<Z5, Column_pairing>,
							Set_column<Z5, Dummy_pairing_5>,
							Unordered_set_column<Z5, Column_pairing>,
							Unordered_set_column<Z5, Dummy_pairing_5>,
							Vector_column<Z5, Column_pairing>,
							Vector_column<Z5, Dummy_pairing_5>
						> list_of_5_columns_without_row;

typedef boost::mpl::list<Z2_heap_column<Column_pairing>,
							Z2_heap_column<Dummy_pairing_2>,
							Z2_list_column<Column_pairing>,
							Z2_list_column<Dummy_pairing_2>,
							Z2_set_column<Column_pairing>,
							Z2_set_column<Dummy_pairing_2>,
							Z2_unordered_set_column<Column_pairing>,
							Z2_unordered_set_column<Dummy_pairing_2>,
							Z2_vector_column<Column_pairing>,
							Z2_vector_column<Dummy_pairing_2>
						> list_of_2_columns_without_row;

typedef boost::mpl::list<List_column<Z5, Column_pairing>,
							Set_column<Z5, Column_pairing>,
							Unordered_set_column<Z5, Column_pairing>,
							Vector_column<Z5, Column_pairing>,
							Z2_heap_column<Column_pairing>,
							Z2_list_column<Column_pairing>,
							Z2_set_column<Column_pairing>,
							Z2_unordered_set_column<Column_pairing>,
							Z2_vector_column<Column_pairing>
						> list_of_columns_with_pairing_without_row;

typedef boost::mpl::list<Reduced_cell_list_column_with_row<List_5_matrix_with_pairing, Z5, Column_pairing>,
							Reduced_cell_set_column_with_row<Set_5_matrix_with_pairing, Z5, Column_pairing>,
							Z2_reduced_cell_list_column_with_row<List_2_matrix_with_pairing, Column_pairing>,
							Z2_reduced_cell_set_column_with_row<Set_2_matrix_with_pairing, Column_pairing>
						> list_of_columns_with_pairing_with_row;

template<class Column>
void common_5_test(std::vector<Column> &matrix)
{
	std::set<std::pair<unsigned int,unsigned int> > rows;
	for (auto c : matrix[0]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	auto it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 0);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 2);
	++it;
	BOOST_CHECK_EQUAL(it->first, 3);
	BOOST_CHECK_EQUAL(it->second, 3);
	++it;
	BOOST_CHECK_EQUAL(it->first, 5);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK(it == rows.end());

	rows.clear();
	for (auto c : matrix[1]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 0);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 2);
	++it;
	BOOST_CHECK_EQUAL(it->first, 2);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 5);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 6);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK(it == rows.end());

	rows.clear();
	for (auto c : matrix[2]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 0);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 3);
	++it;
	BOOST_CHECK_EQUAL(it->first, 2);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 5);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 6);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 5);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 6);
	BOOST_CHECK_EQUAL(matrix[2].get_pivot(), 6);
	BOOST_CHECK_EQUAL(matrix[0].get_pivot_value(), 4u);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot_value(), 1u);
	BOOST_CHECK_EQUAL(matrix[2].get_pivot_value(), 4u);
	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[2].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(!matrix[1].is_empty());
	BOOST_CHECK(!matrix[2].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(1));
	BOOST_CHECK(matrix[1].is_non_zero(2));
	BOOST_CHECK(!matrix[2].is_non_zero(3));

	matrix[0] += matrix[1];
	matrix[1] += matrix[2];

	rows.clear();
	for (auto c : matrix[0]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 2);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 3);
	BOOST_CHECK_EQUAL(it->second, 3);
	++it;
	BOOST_CHECK_EQUAL(it->first, 6);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK(matrix[1].begin() == matrix[1].end());

//	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 5);
//	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 6);
//	BOOST_CHECK_EQUAL(matrix[0].get_pivot_value(), 4u);
//	BOOST_CHECK_EQUAL(matrix[1].get_pivot_value(), 1u);
	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(matrix[1].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(1));
	BOOST_CHECK(!matrix[1].is_non_zero(2));

	matrix[0] *= 4;
	matrix[1] *= 2;
	matrix[2] *= 0;

	rows.clear();
	for (auto c : matrix[0]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 2);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 3);
	BOOST_CHECK_EQUAL(it->second, 2);
	++it;
	BOOST_CHECK_EQUAL(it->first, 6);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK(matrix[1].begin() == matrix[1].end());
	BOOST_CHECK(matrix[2].begin() == matrix[2].end());

//	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 5);
//	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 6);
//	BOOST_CHECK_EQUAL(matrix[0].get_pivot_value(), 4u);
//	BOOST_CHECK_EQUAL(matrix[1].get_pivot_value(), 1u);
	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(matrix[1].is_empty());
	BOOST_CHECK(matrix[2].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(1));
	BOOST_CHECK(!matrix[1].is_non_zero(2));
	BOOST_CHECK(!matrix[2].is_non_zero(3));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_common_5_without_row, Column, list_of_5_columns_without_row) {
	std::vector<Column> matrix;

	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}, 4));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(4)},{1,Z5(2)},{2,Z5(1)},{5,Z5(1)},{6,Z5(1)}}, 4));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(3)},{2,Z5(4)},{5,Z5(4)},{6,Z5(4)}}, 4));

	common_5_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_common_5_with_row, Column, list_of_5_columns_with_row) {
	std::vector<Column> matrix;
	std::vector<unsigned int> pivotToColumnIndex{0,0,0,0,0,0,1,2};

	matrix.push_back(Column(0, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}, 4, matrix, pivotToColumnIndex));
	matrix.push_back(Column(1, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(4)},{1,Z5(2)},{2,Z5(1)},{5,Z5(1)},{6,Z5(1)}}, 4, matrix, pivotToColumnIndex));
	matrix.push_back(Column(2, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(3)},{2,Z5(4)},{5,Z5(4)},{6,Z5(4)}}, 4, matrix, pivotToColumnIndex));

	common_5_test(matrix);
}

template<class Column>
void common_2_test(std::vector<Column> &matrix)
{
	std::set<unsigned int> rows;
	for (auto c : matrix[0]){
		rows.insert(c.get_row_index());
	}
	auto it = rows.begin();
	BOOST_CHECK_EQUAL(*it, 0);
	++it;
	BOOST_CHECK_EQUAL(*it, 1);
	++it;
	BOOST_CHECK_EQUAL(*it, 3);
	++it;
	BOOST_CHECK_EQUAL(*it, 5);
	++it;
	BOOST_CHECK(it == rows.end());

	rows.clear();
	for (auto c : matrix[1]){
		rows.insert(c.get_row_index());
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(*it, 0);
	++it;
	BOOST_CHECK_EQUAL(*it, 1);
	++it;
	BOOST_CHECK_EQUAL(*it, 2);
	++it;
	BOOST_CHECK_EQUAL(*it, 5);
	++it;
	BOOST_CHECK_EQUAL(*it, 6);
	++it;
	BOOST_CHECK(it == rows.end());

	rows.clear();
	for (auto c : matrix[2]){
		rows.insert(c.get_row_index());
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(*it, 0);
	++it;
	BOOST_CHECK_EQUAL(*it, 1);
	++it;
	BOOST_CHECK_EQUAL(*it, 2);
	++it;
	BOOST_CHECK_EQUAL(*it, 5);
	++it;
	BOOST_CHECK_EQUAL(*it, 6);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 5);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 6);
	BOOST_CHECK_EQUAL(matrix[2].get_pivot(), 6);
	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[2].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(!matrix[1].is_empty());
	BOOST_CHECK(!matrix[2].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(1));
	BOOST_CHECK(matrix[1].is_non_zero(2));
	BOOST_CHECK(!matrix[2].is_non_zero(3));

	matrix[0] += matrix[1];
	matrix[1] += matrix[2];

	rows.clear();
	for (auto c : matrix[0]){
		rows.insert(c.get_row_index());
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(*it, 2);
	++it;
	BOOST_CHECK_EQUAL(*it, 3);
	++it;
	BOOST_CHECK_EQUAL(*it, 6);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK(matrix[1].begin() == matrix[1].end());

//	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 5);
//	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 6);
	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(matrix[1].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(2));
	BOOST_CHECK(!matrix[1].is_non_zero(2));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_common_2_without_row, Column, list_of_2_columns_without_row) {
	std::vector<Column> matrix;

	matrix.push_back(Column(std::vector<unsigned int>{0,1,3,5}, 4));
	matrix.push_back(Column(std::vector<unsigned int>{0,1,2,5,6}, 4));
	matrix.push_back(Column(std::vector<unsigned int>{0,1,2,5,6}, 4));

	common_2_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_common_2_with_row, Column, list_of_2_columns_with_row) {
	std::vector<Column> matrix;
	std::vector<unsigned int> pivotToColumnIndex{0,0,0,0,0,0,1,2};

	matrix.push_back(Column(0, std::vector<unsigned int>{0,1,3,5}, 4, matrix, pivotToColumnIndex));
	matrix.push_back(Column(1, std::vector<unsigned int>{0,1,2,5,6}, 4, matrix, pivotToColumnIndex));
	matrix.push_back(Column(2, std::vector<unsigned int>{0,1,2,5,6}, 4, matrix, pivotToColumnIndex));

	common_2_test(matrix);
}

template<class Column>
void pairing_test(std::vector<Column> &matrix)
{
	matrix[0].assign_paired_chain(1);
	matrix[1].assign_paired_chain(0);
	matrix[3].assign_paired_chain(2);

	BOOST_CHECK(matrix[0].is_paired());
	BOOST_CHECK(matrix[1].is_paired());
	BOOST_CHECK(!matrix[2].is_paired());
	BOOST_CHECK(matrix[3].is_paired());
	BOOST_CHECK_EQUAL(matrix[0].get_paired_chain_index(), 1);
	BOOST_CHECK_EQUAL(matrix[1].get_paired_chain_index(), 0);
	BOOST_CHECK_EQUAL(matrix[2].get_paired_chain_index(), -1);
	BOOST_CHECK_EQUAL(matrix[3].get_paired_chain_index(), 2);

	matrix[3].unassign_paired_chain();
	matrix[1].unassign_paired_chain();

	BOOST_CHECK(matrix[0].is_paired());
	BOOST_CHECK(!matrix[1].is_paired());
	BOOST_CHECK(!matrix[2].is_paired());
	BOOST_CHECK(!matrix[3].is_paired());
	BOOST_CHECK_EQUAL(matrix[0].get_paired_chain_index(), 1);
	BOOST_CHECK_EQUAL(matrix[1].get_paired_chain_index(), -1);
	BOOST_CHECK_EQUAL(matrix[2].get_paired_chain_index(), -1);
	BOOST_CHECK_EQUAL(matrix[3].get_paired_chain_index(), -1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_pairing_option_without_row, Column, list_of_columns_with_pairing_without_row) {
	std::vector<Column> matrix(4);
	pairing_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_pairing_option_with_row, Column, list_of_columns_with_pairing_with_row) {
	std::vector<Column> matrix;
	std::vector<unsigned int> pivotToColumnIndex{0,0,0,0};
	for (int i = 0; i < 4; ++i)
		matrix.push_back(Column(matrix, pivotToColumnIndex));

	pairing_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_without_row_methods, Column, list_of_5_columns_without_row) {
	Column col1;
	Column col2(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}});
	Column col3(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(4)},{1,Z5(2)},{2,Z5(1)},{5,Z5(1)},{6,Z5(1)}}, 4);
	Column col4(col3);

	BOOST_CHECK_EQUAL(col2.get_dimension(), 3);
	BOOST_CHECK_EQUAL(col3.get_dimension(), 4);

	auto res = col1.get_content(7);
	for (auto& f : res) BOOST_CHECK_EQUAL(f, 0u);

	res = col2.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 0u);
	BOOST_CHECK_EQUAL(res[3], 3u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 4u);
	BOOST_CHECK_EQUAL(res[6], 0u);

	res = col3.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 4u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 0u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 1u);
	BOOST_CHECK_EQUAL(res[6], 1u);

	BOOST_CHECK(res == col4.get_content(7));

	col4.clear(5);
	res = col4.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 4u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 0u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 0u);
	BOOST_CHECK_EQUAL(res[6], 1u);

	col4.clear();
	BOOST_CHECK(col4.is_empty());
	BOOST_CHECK(col1.get_content(7) == col4.get_content(7));

	std::vector<unsigned int> permutation{0,5,2,3,1,4,6};
	col3.reorder(permutation);
	res = col3.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 4u);
	BOOST_CHECK_EQUAL(res[1], 0u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 0u);
	BOOST_CHECK_EQUAL(res[4], 1u);
	BOOST_CHECK_EQUAL(res[5], 2u);
	BOOST_CHECK_EQUAL(res[6], 1u);

	Column add = col2 + col3;
	res = add.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 0u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 3u);
	BOOST_CHECK_EQUAL(res[4], 1u);
	BOOST_CHECK_EQUAL(res[5], 1u);
	BOOST_CHECK_EQUAL(res[6], 1u);

	Column mul1 = col2 * 2;
	res = mul1.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 2u);
	BOOST_CHECK_EQUAL(res[1], 4u);
	BOOST_CHECK_EQUAL(res[2], 0u);
	BOOST_CHECK_EQUAL(res[3], 1u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 3u);
	BOOST_CHECK_EQUAL(res[6], 0u);

	Column mul2 = 2 * col2;
	BOOST_CHECK(res == mul2.get_content(7));

	swap(col1, col2);
	res = col2.get_content(7);
	for (auto& f : res) BOOST_CHECK_EQUAL(f, 0u);
	res = col1.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 0u);
	BOOST_CHECK_EQUAL(res[3], 3u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 4u);
	BOOST_CHECK_EQUAL(res[6], 0u);

	Column move(std::move(col3));
	res = move.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 4u);
	BOOST_CHECK_EQUAL(res[1], 0u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 0u);
	BOOST_CHECK_EQUAL(res[4], 1u);
	BOOST_CHECK_EQUAL(res[5], 2u);
	BOOST_CHECK_EQUAL(res[6], 1u);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_z2_without_row_methods, Column, list_of_2_columns_without_row) {
	Column col1;
	Column col2(std::vector<unsigned int>{0,1,3,5});
	Column col3(std::vector<unsigned int>{0,1,2,5,6}, 4);
	Column col4(col3);

	BOOST_CHECK_EQUAL(col2.get_dimension(), 3);
	BOOST_CHECK_EQUAL(col3.get_dimension(), 4);

	auto res = col1.get_content(7);
	for (auto f : res) BOOST_CHECK_EQUAL(f, 0);

	res = col2.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 0);
	BOOST_CHECK_EQUAL(res[3], 1);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 0);

	res = col3.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 0);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 1);

	BOOST_CHECK(res == col4.get_content(7));

	col4.clear(5);
	res = col4.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 0);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 0);
	BOOST_CHECK_EQUAL(res[6], 1);

	col4.clear();
	BOOST_CHECK(col4.is_empty());
	BOOST_CHECK(col1.get_content(7) == col4.get_content(7));

	std::vector<unsigned int> permutation{0,5,2,3,1,4,6};
	col3.reorder(permutation);
	res = col3.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 0);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 0);
	BOOST_CHECK_EQUAL(res[4], 1);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 1);

	Column add = col2 + col3;
	res = add.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 0);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 1);
	BOOST_CHECK_EQUAL(res[4], 1);
	BOOST_CHECK_EQUAL(res[5], 0);
	BOOST_CHECK_EQUAL(res[6], 1);

	swap(col1, col2);
	res = col2.get_content(7);
	for (auto f : res) BOOST_CHECK_EQUAL(f, 0u);
	res = col1.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 0);
	BOOST_CHECK_EQUAL(res[3], 1);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 0);

	Column move(std::move(col3));
	res = move.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 0);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 0);
	BOOST_CHECK_EQUAL(res[4], 1);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_methods, Column, list_of_5_columns_with_row) {
	std::vector<Column> matrix;
	std::vector<unsigned int> pivotToColumnIndex{0,1,2,3};

	std::vector<std::vector<std::pair<unsigned int,Z5> > > chains;
	chains.push_back({{0,Z5(1)}});
	chains.push_back({{1,Z5(3)}});
	chains.push_back({{1,Z5(2)},{2,Z5(4)}});
	chains.push_back({{0,Z5(2)},{2,Z5(3)},{3,Z5(1)}});
	chains.emplace_back();
	chains.push_back({{0,Z5(1)}});

	std::vector<std::vector<std::pair<unsigned int,Z5> > > rows;
	rows.push_back({{0,Z5(1)},{3,Z5(2)}});
	rows.push_back({{1,Z5(3)},{2,Z5(2)}});
	rows.push_back({{2,Z5(4)},{3,Z5(3)}});
	rows.push_back({{3,Z5(1)}});
	rows.emplace_back();
	rows.push_back({{0,Z5(1)},{3,Z5(2)}});

	for (unsigned int i = 0; i < 4; ++i){
		matrix.push_back(Column(i, chains[i], 4, matrix, pivotToColumnIndex));
		for (auto& cell : matrix.back().get_column()){
			matrix[cell.get_row_index()].get_row().push_back(cell);
		}
	}
	matrix.push_back(Column(matrix, pivotToColumnIndex));
	matrix.push_back(matrix[0]);

	for (unsigned int i = 0; i < matrix.size(); ++i){
		auto& col = matrix[i].get_column();
		auto itCol = col.begin();
		auto itChain = chains[i].begin();
		while (itChain != chains[i].end()){
			BOOST_CHECK_EQUAL(itCol->get_row_index(), itChain->first);
			BOOST_CHECK_EQUAL(itCol->get_element(), itChain->second);
			++itCol; ++itChain;
		}

		auto& row = matrix[i].get_row();
		auto itRow = row.begin();
		itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(itRow->get_column_index(), itChain->first);
			BOOST_CHECK_EQUAL(itRow->get_element(), itChain->second);
			++itRow; ++itChain;
		}
	}

	BOOST_CHECK_EQUAL(matrix[0].get_lowest_simplex_index(), 0);
	BOOST_CHECK_EQUAL(matrix[1].get_lowest_simplex_index(), 1);
	BOOST_CHECK_EQUAL(matrix[2].get_lowest_simplex_index(), 2);
	BOOST_CHECK_EQUAL(matrix[3].get_lowest_simplex_index(), 3);
	BOOST_CHECK_EQUAL(matrix[4].get_lowest_simplex_index(), -1);
	BOOST_CHECK_EQUAL(matrix[5].get_lowest_simplex_index(), 0);

	matrix[1].swap_rows(matrix[2]);
	auto itRow = matrix[1].get_row().begin();
	auto itChain = rows[2].begin();
	while (itChain != rows[2].end()){
		BOOST_CHECK_EQUAL(itRow->get_column_index(), itChain->first);
		BOOST_CHECK_EQUAL(itRow->get_element(), itChain->second);
		++itRow; ++itChain;
	}
	itRow = matrix[2].get_row().begin();
	itChain = rows[1].begin();
	while (itChain != rows[1].end()){
		BOOST_CHECK_EQUAL(itRow->get_column_index(), itChain->first);
		BOOST_CHECK_EQUAL(itRow->get_element(), itChain->second);
		++itRow; ++itChain;
	}

	matrix[1].swap_lowest_simplex_index(matrix[2]);
	BOOST_CHECK_EQUAL(matrix[1].get_lowest_simplex_index(), 2);
	BOOST_CHECK_EQUAL(matrix[2].get_lowest_simplex_index(), 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_z2_with_row_methods, Column, list_of_2_columns_with_row) {
	std::vector<Column> matrix;
	std::vector<unsigned int> pivotToColumnIndex{0,1,2,3};

	std::vector<std::vector<unsigned int> > chains;
	chains.push_back({0});
	chains.push_back({1});
	chains.push_back({1,2});
	chains.push_back({0,2,3});
	chains.emplace_back();
	chains.push_back({0});

	std::vector<std::vector<unsigned int> > rows;
	rows.push_back({0,3});
	rows.push_back({1,2});
	rows.push_back({2,3});
	rows.push_back({3});
	rows.emplace_back();
	rows.push_back({0,3});

	for (unsigned int i = 0; i < 4; ++i){
		matrix.push_back(Column(i, chains[i], 4, matrix, pivotToColumnIndex));
		for (auto& cell : matrix.back().get_column()){
			matrix[cell.get_row_index()].get_row().push_back(cell);
		}
	}
	matrix.push_back(Column(matrix, pivotToColumnIndex));
	matrix.push_back(matrix[0]);

	for (unsigned int i = 0; i < matrix.size(); ++i){
		auto& col = matrix[i].get_column();
		auto itCol = col.begin();
		auto itChain = chains[i].begin();
		while (itChain != chains[i].end()){
			BOOST_CHECK_EQUAL(itCol->get_row_index(), *itChain);
			++itCol; ++itChain;
		}

		auto& row = matrix[i].get_row();
		auto itRow = row.begin();
		itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(itRow->get_column_index(), *itChain);
			++itRow; ++itChain;
		}
	}

	BOOST_CHECK_EQUAL(matrix[0].get_lowest_simplex_index(), 0);
	BOOST_CHECK_EQUAL(matrix[1].get_lowest_simplex_index(), 1);
	BOOST_CHECK_EQUAL(matrix[2].get_lowest_simplex_index(), 2);
	BOOST_CHECK_EQUAL(matrix[3].get_lowest_simplex_index(), 3);
	BOOST_CHECK_EQUAL(matrix[4].get_lowest_simplex_index(), -1);
	BOOST_CHECK_EQUAL(matrix[5].get_lowest_simplex_index(), 0);

	matrix[1].swap_rows(matrix[2]);
	auto itRow = matrix[1].get_row().begin();
	auto itChain = rows[2].begin();
	while (itChain != rows[2].end()){
		BOOST_CHECK_EQUAL(itRow->get_column_index(), *itChain);
		++itRow; ++itChain;
	}
	itRow = matrix[2].get_row().begin();
	itChain = rows[1].begin();
	while (itChain != rows[1].end()){
		BOOST_CHECK_EQUAL(itRow->get_column_index(), *itChain);
		++itRow; ++itChain;
	}

	matrix[1].swap_lowest_simplex_index(matrix[2]);
	BOOST_CHECK_EQUAL(matrix[1].get_lowest_simplex_index(), 2);
	BOOST_CHECK_EQUAL(matrix[2].get_lowest_simplex_index(), 1);
}


