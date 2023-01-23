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
#define BOOST_MPL_LIMIT_LIST_SIZE 50
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>
#include <random>
#include <vector>
#include <set>
#include <utility>

#include "gudhi/matrix.h"
#include "gudhi/options.h"
#include "gudhi/utilities/Z2_field.h"
#include "gudhi/utilities/Zp_field.h"
#include "gudhi/utilities/utilities.h"

using Gudhi::persistence_matrix::Z2_field_element;
using Gudhi::persistence_matrix::Zp_field_element;
using Gudhi::persistence_matrix::Matrix;
using Gudhi::persistence_matrix::Representative_cycles_options;
using Gudhi::persistence_matrix::Default_options;
using Gudhi::persistence_matrix::Zigzag_options;
using Gudhi::persistence_matrix::Multi_persistence_options;
using Gudhi::persistence_matrix::Cohomology_persistence_options;
using Gudhi::persistence_matrix::Column_types;

using Z5 = Zp_field_element<5>;
using Z2 = Zp_field_element<2>;

template<class Field_type, Column_types column_type>
struct opt_ra_i_r : Default_options<Field_type, column_type, false>{
	static const bool has_column_compression = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_r : Default_options<Field_type, column_type, false>{
	static const bool has_column_compression = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra : Default_options<Field_type, column_type, false>{
	static const bool has_column_compression = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_i : Default_options<Field_type, column_type, false>{
	static const bool has_column_compression = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = false;
};

template<class Field_type, Column_types column_type>
struct opt_r : Default_options<Field_type, column_type, false>{
	static const bool has_column_compression = false;
	static const bool has_row_access = false;
	static const bool has_removable_columns = true;
};

template<class Field_type, Column_types column_type>
struct opt : Default_options<Field_type, column_type, false>{
	static const bool has_column_compression = false;
	static const bool has_row_access = false;
	static const bool has_removable_columns = false;
};

template<class Field_type, Column_types column_type>
struct opt_cc_ra : Default_options<Field_type, column_type, false>{
	static const bool has_column_compression = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_cc_ra_i : Default_options<Field_type, column_type, false>{
	static const bool has_column_compression = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
};

template<class Field_type, Column_types column_type>
struct opt_cc : Default_options<Field_type, column_type, false>{
	static const bool has_column_compression = true;
	static const bool has_row_access = false;
};

void build_boundary_matrix(std::vector<std::vector<unsigned int> >& boundaries)
{
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.push_back(std::vector<unsigned int>{0,1,4});
	boundaries.emplace_back();
	boundaries.push_back(std::vector<unsigned int>{1,2});
	boundaries.push_back(std::vector<unsigned int>{0,4});
	boundaries.push_back(std::vector<unsigned int>{3,4,5});
	boundaries.push_back(std::vector<unsigned int>{1,2});
	boundaries.push_back(std::vector<unsigned int>{0,1,4});
	boundaries.push_back(std::vector<unsigned int>{0,4});
	boundaries.push_back(std::vector<unsigned int>{0,1,4});
}

template<typename Field_type>
void build_boundary_matrix(std::vector<std::vector<std::pair<unsigned int,Field_type> > >& boundaries)
{
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{0,Field_type(1)},{1,Field_type(4)},{4,Field_type(1)}});
	boundaries.emplace_back();
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{1,Field_type(1)},{2,Field_type(4)}});
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{0,Field_type(1)},{4,Field_type(4)}});
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{3,Field_type(1)},{4,Field_type(1)},{5,Field_type(4)}});
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{1,Field_type(1)},{2,Field_type(4)}});
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{0,Field_type(1)},{1,Field_type(4)},{4,Field_type(1)}});
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{0,Field_type(1)},{4,Field_type(4)}});
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{0,Field_type(1)},{1,Field_type(4)},{4,Field_type(1)}});
}

typedef boost::mpl::list<Matrix<opt_ra_i_r<Z5,Column_types::LIST> >,
							Matrix<opt_ra_i_r<Z2,Column_types::LIST> >,
							Matrix<opt_ra_i_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_i_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::LIST> >,
							Matrix<opt_ra_r<Z2,Column_types::LIST> >,
							Matrix<opt_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_r<Z2,Column_types::SET> >,
							Matrix<opt_ra_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z5,Column_types::LIST> >,
							Matrix<opt_r<Z2,Column_types::LIST> >,
							Matrix<opt_r<Z5,Column_types::SET> >,
							Matrix<opt_r<Z2,Column_types::SET> >,
							Matrix<opt_r<Z2,Column_types::HEAP> >,
							Matrix<opt_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_matrix_types_with_removable_columns;

typedef boost::mpl::list<Matrix<opt_ra_i<Z5,Column_types::LIST> >,
							Matrix<opt_ra_i<Z2,Column_types::LIST> >,
							Matrix<opt_ra_i<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_i<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z5,Column_types::LIST> >,
							Matrix<opt_ra<Z2,Column_types::LIST> >,
							Matrix<opt_ra<Z5,Column_types::SET> >,
							Matrix<opt_ra<Z2,Column_types::SET> >,
							Matrix<opt_ra<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z5,Column_types::LIST> >,
							Matrix<opt<Z2,Column_types::LIST> >,
							Matrix<opt<Z5,Column_types::SET> >,
							Matrix<opt<Z2,Column_types::SET> >,
							Matrix<opt<Z2,Column_types::HEAP> >,
							Matrix<opt<Z5,Column_types::VECTOR> >,
							Matrix<opt<Z2,Column_types::VECTOR> >,
							Matrix<opt<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_matrix_types;

typedef boost::mpl::list<Matrix<opt_cc_ra_i<Z5,Column_types::LIST> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::LIST> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::VECTOR> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::VECTOR> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc_ra<Z5,Column_types::LIST> >,
							Matrix<opt_cc_ra<Z2,Column_types::LIST> >,
							Matrix<opt_cc_ra<Z5,Column_types::SET> >,
							Matrix<opt_cc_ra<Z2,Column_types::SET> >,
							Matrix<opt_cc_ra<Z5,Column_types::VECTOR> >,
							Matrix<opt_cc_ra<Z2,Column_types::VECTOR> >,
							Matrix<opt_cc_ra<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc_ra<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc_ra<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc_ra<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc<Z5,Column_types::LIST> >,
							Matrix<opt_cc<Z2,Column_types::LIST> >,
							Matrix<opt_cc<Z5,Column_types::SET> >,
							Matrix<opt_cc<Z2,Column_types::SET> >,
							Matrix<opt_cc<Z5,Column_types::VECTOR> >,
							Matrix<opt_cc<Z2,Column_types::VECTOR> >,
							Matrix<opt_cc<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_matrix_types_with_column_compression;

template<class Matrix>
void test_constructors(){
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	//default constructor
	Matrix m;
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 0);

	//constructor from given matrix
	Matrix mb(ordered_boundaries);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 12);

	//constructor reserving column space
	Matrix mr(5);
	BOOST_CHECK_EQUAL(mr.get_number_of_columns(), 0);

	//copy constructor
	Matrix mc1(mb);
	Matrix mc2 = mb;
	BOOST_CHECK_EQUAL(mc1.get_number_of_columns(), 12);
	BOOST_CHECK_EQUAL(mc2.get_number_of_columns(), 12);

	//move constructor
	Matrix mm(std::move(mb));
	BOOST_CHECK_EQUAL(mm.get_number_of_columns(), 12);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 0);

	//swap
	swap(mm, mb);
	BOOST_CHECK_EQUAL(mm.get_number_of_columns(), 0);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 12);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_constructors_removable_columns, Matrix, list_of_matrix_types_with_removable_columns) {
	test_constructors<Matrix>();
}
BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_constructors, Matrix, list_of_matrix_types) {
	test_constructors<Matrix>();
}
BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_constructors_with_column_compression, Matrix, list_of_matrix_types_with_column_compression) {
	test_constructors<Matrix>();
}

template<class Matrix>
void test_basic_methods(){
	using boundary_type = typename Matrix::boundary_type;
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	boundary_type boundary2 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	boundary_type boundary1 = ordered_boundaries.back();
	ordered_boundaries.pop_back();

	Matrix m(ordered_boundaries);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(m.is_zero_column(4));
	BOOST_CHECK(!m.is_zero_column(5));
	BOOST_CHECK(!m.is_zero_column(6));
	BOOST_CHECK(!m.is_zero_column(7));
	BOOST_CHECK(!m.is_zero_column(8));
	BOOST_CHECK(!m.is_zero_column(9));

	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 10);
	m.insert_boundary(boundary1);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 11);
	m.insert_boundary(boundary2);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 12);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(m.is_zero_column(4));
	BOOST_CHECK(!m.is_zero_column(5));
	BOOST_CHECK(!m.is_zero_column(6));
	BOOST_CHECK(!m.is_zero_column(7));
	BOOST_CHECK(!m.is_zero_column(8));
	BOOST_CHECK(!m.is_zero_column(9));
	BOOST_CHECK(!m.is_zero_column(10));
	BOOST_CHECK(!m.is_zero_column(11));
	BOOST_CHECK(!m.is_zero_cell(3,1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_basics_removable_columns, Matrix, list_of_matrix_types_with_removable_columns) {
	test_basic_methods<Matrix>();
}
BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_basics, Matrix, list_of_matrix_types) {
	test_basic_methods<Matrix>();
}
BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_basics_with_column_compression, Matrix, list_of_matrix_types_with_column_compression) {
	test_basic_methods<Matrix>();
}

typedef boost::mpl::list<Matrix<opt_ra_i_r<Z2,Column_types::LIST> >,
							Matrix<opt_ra_i_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z2,Column_types::LIST> >,
							Matrix<opt_ra_r<Z2,Column_types::SET> >,
							Matrix<opt_ra_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z2,Column_types::LIST> >,
							Matrix<opt_r<Z2,Column_types::SET> >,
							Matrix<opt_r<Z2,Column_types::HEAP> >,
							Matrix<opt_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i<Z2,Column_types::LIST> >,
							Matrix<opt_ra_i<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z2,Column_types::LIST> >,
							Matrix<opt_ra<Z2,Column_types::SET> >,
							Matrix<opt_ra<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z2,Column_types::LIST> >,
							Matrix<opt<Z2,Column_types::SET> >,
							Matrix<opt<Z2,Column_types::HEAP> >,
							Matrix<opt<Z2,Column_types::VECTOR> >,
							Matrix<opt<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_z2_matrix_types_without_column_compression;

typedef boost::mpl::list<Matrix<opt_cc_ra_i<Z2,Column_types::LIST> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::VECTOR> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc_ra<Z2,Column_types::LIST> >,
							Matrix<opt_cc_ra<Z2,Column_types::SET> >,
							Matrix<opt_cc_ra<Z2,Column_types::VECTOR> >,
							Matrix<opt_cc_ra<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc_ra<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc<Z2,Column_types::LIST> >,
							Matrix<opt_cc<Z2,Column_types::SET> >,
							Matrix<opt_cc<Z2,Column_types::VECTOR> >,
							Matrix<opt_cc<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_z2_matrix_types_with_column_compression;

template<class Column_type>
void test_equal_columns(const std::vector<unsigned int>& col1, const Column_type& col2)
{
	std::set<unsigned int> rowIndices;
	//not all columns are ordered
	for (const auto& cell : col2){
		rowIndices.insert(cell.get_row_index());
	}
	auto itCol = rowIndices.begin();
	auto itOrdB = col1.begin();
	while (itCol != rowIndices.end()) {
		BOOST_CHECK_EQUAL(*itCol, *itOrdB);
		itCol++; itOrdB++;
	}
}

template<class Column_type>
void test_equal_columns(const std::vector<std::pair<unsigned int,Z5> >& col1, const Column_type& col2)
{
	std::set<std::pair<unsigned int,unsigned int> > rowIndices;
	//not all columns are ordered
	for (const auto& cell : col2){
		rowIndices.insert({cell.get_row_index(), cell.get_element()});
	}
	auto itCol = rowIndices.begin();
	auto itOrdB = col1.begin();
	while (itCol != rowIndices.end()) {
		BOOST_CHECK_EQUAL(itCol->first, itOrdB->first);
		BOOST_CHECK_EQUAL(itCol->second, itOrdB->second);
		itCol++; itOrdB++;
	}
}

template<class Matrix>
void test_z2_insertion_and_addition(){
	using boundary_type = typename Matrix::boundary_type;

	std::vector<std::vector<unsigned int> > ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	boundary_type boundary2 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	boundary_type boundary1 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	Matrix m(ordered_boundaries);

	unsigned int i = 0;
	for (auto& b : ordered_boundaries){
		const auto& col = m.get_column(i++);
		test_equal_columns(b, col);
	}

	m.insert_boundary(boundary1);
	ordered_boundaries.push_back(boundary1);
	m.insert_boundary(boundary2);
	ordered_boundaries.push_back(boundary2);

	i = 0;
	for (auto& b : ordered_boundaries){
		const auto& col = m.get_column(i++);
		test_equal_columns(b, col);
	}

	m.add_to(6, 7);

	std::vector<bool> column = m.get_column(7).get_content();
	BOOST_CHECK_EQUAL(column.size(), 6);
	auto it = column.begin();
	BOOST_CHECK_EQUAL(*it++, 1);
	BOOST_CHECK_EQUAL(*it++, 0);
	BOOST_CHECK_EQUAL(*it++, 0);
	BOOST_CHECK_EQUAL(*it++, 1);
	BOOST_CHECK_EQUAL(*it++, 0);
	BOOST_CHECK_EQUAL(*it, 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_z2_insertion_and_addition, Matrix, list_of_z2_matrix_types_without_column_compression) {
	test_z2_insertion_and_addition<Matrix>();
}
BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_z2_insertion_and_addition_with_cc, Matrix, list_of_z2_matrix_types_with_column_compression) {
	test_z2_insertion_and_addition<Matrix>();
}

typedef boost::mpl::list<Matrix<opt_ra_i_r<Z5,Column_types::LIST> >,
							Matrix<opt_ra_i_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::LIST> >,
							Matrix<opt_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z5,Column_types::LIST> >,
							Matrix<opt_r<Z5,Column_types::SET> >,
							Matrix<opt_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i<Z5,Column_types::LIST> >,
							Matrix<opt_ra_i<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z5,Column_types::LIST> >,
							Matrix<opt_ra<Z5,Column_types::SET> >,
							Matrix<opt_ra<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z5,Column_types::LIST> >,
							Matrix<opt<Z5,Column_types::SET> >,
							Matrix<opt<Z5,Column_types::VECTOR> >,
							Matrix<opt<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt<Z5,Column_types::INTRUSIVE_SET> >
						> list_of_z5_matrix_types_without_column_compression;

typedef boost::mpl::list<Matrix<opt_cc_ra_i<Z5,Column_types::LIST> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::VECTOR> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc_ra<Z5,Column_types::LIST> >,
							Matrix<opt_cc_ra<Z5,Column_types::SET> >,
							Matrix<opt_cc_ra<Z5,Column_types::VECTOR> >,
							Matrix<opt_cc_ra<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc_ra<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc<Z5,Column_types::LIST> >,
							Matrix<opt_cc<Z5,Column_types::SET> >,
							Matrix<opt_cc<Z5,Column_types::VECTOR> >,
							Matrix<opt_cc<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc<Z5,Column_types::INTRUSIVE_SET> >
						> list_of_z5_matrix_types_with_column_compression;

template<class Matrix>
void test_z5_insertion_and_addition(){
	using boundary_type = typename Matrix::boundary_type;

	std::vector<std::vector<std::pair<unsigned int,Z5> > > ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	boundary_type boundary2 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	boundary_type boundary1 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	Matrix m(ordered_boundaries);

	unsigned int i = 0;
	for (auto& b : ordered_boundaries){
		const auto& col = m.get_column(i++);
		test_equal_columns(b, col);
	}

	m.insert_boundary(boundary1);
	ordered_boundaries.push_back(boundary1);
	m.insert_boundary(boundary2);
	ordered_boundaries.push_back(boundary2);

	i = 0;
	for (auto& b : ordered_boundaries){
		const auto& col = m.get_column(i++);
		test_equal_columns(b, col);
	}

	m.add_to(6, 7);

	std::set<std::pair<unsigned int,unsigned int> > rowIndices;
	//not all columns are ordered
	for (auto& cell : m.get_column(7)){
		rowIndices.insert({cell.get_row_index(), cell.get_element()});
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 3);
	auto it = rowIndices.begin();
	BOOST_CHECK_EQUAL(it->first, 0);
	BOOST_CHECK_EQUAL((it++)->second, 1);
	BOOST_CHECK_EQUAL(it->first, 3);
	BOOST_CHECK_EQUAL((it++)->second, 1);
	BOOST_CHECK_EQUAL(it->first, 5);
	BOOST_CHECK_EQUAL(it->second, 4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_z5_insertion_and_addition, Matrix, list_of_z5_matrix_types_without_column_compression) {
	test_z5_insertion_and_addition<Matrix>();
}
BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_z5_insertion_and_addition_with_cc, Matrix, list_of_z5_matrix_types_with_column_compression) {
	test_z5_insertion_and_addition<Matrix>();
}

typedef boost::mpl::list<Matrix<opt_ra_i_r<Z5,Column_types::LIST> >,
							Matrix<opt_ra_i_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::LIST> >,
							Matrix<opt_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i<Z5,Column_types::LIST> >,
							Matrix<opt_ra_i<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z5,Column_types::LIST> >,
							Matrix<opt_ra<Z5,Column_types::SET> >,
							Matrix<opt_ra<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::LIST> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::VECTOR> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc_ra<Z5,Column_types::LIST> >,
							Matrix<opt_cc_ra<Z5,Column_types::SET> >,
							Matrix<opt_cc_ra<Z5,Column_types::VECTOR> >,
							Matrix<opt_cc_ra<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc_ra<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra<Z5,Column_types::INTRUSIVE_SET> >
						> list_of_z5_matrix_types_with_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_row_access_columns_options, Matrix, list_of_z5_matrix_types_with_row_access) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	std::vector<std::vector<std::pair<unsigned int,Z5> > > rows;
	if constexpr (Matrix::Option_list::has_column_compression){
		//if the union find structure changes, the column_index values of de cells could also change. Change the test with all possibilities?
		rows.push_back({{9,Z5(1)},{10,Z5(1)}});
		rows.push_back({{8,Z5(1)},{9,Z5(4)}});
		rows.push_back({{8,Z5(4)}});
		rows.push_back({{7,Z5(1)}});
		rows.push_back({{7,Z5(1)},{9,Z5(1)},{10,Z5(4)}});
		rows.push_back({{7,Z5(4)}});
	} else {
		rows.push_back({{3,Z5(1)},{6,Z5(1)},{9,Z5(1)},{10,Z5(1)},{11,Z5(1)}});
		rows.push_back({{3,Z5(4)},{5,Z5(1)},{8,Z5(1)},{9,Z5(4)},{11,Z5(4)}});
		rows.push_back({{5,Z5(4)},{8,Z5(4)}});
		rows.push_back({{7,Z5(1)}});
		rows.push_back({{3,Z5(1)},{6,Z5(4)},{7,Z5(1)},{9,Z5(1)},{10,Z5(4)},{11,Z5(1)}});
		rows.push_back({{7,Z5(4)}});
	}

	//rows are unordered
	std::vector<std::set<std::pair<unsigned int,Z5> > > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		for (const auto& cell : mb.get_row(i)){
			ordered_rows[i].insert({cell.get_column_index(), cell.get_element()});
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		auto& row = ordered_rows[i];
		BOOST_CHECK_EQUAL(row.size(), rows[i].size());
		auto itRow = row.begin();
		auto itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(itRow->first, itChain->first);
			BOOST_CHECK_EQUAL(itRow->second, itChain->second);
			++itRow; ++itChain;
		}
	}
}

typedef boost::mpl::list<Matrix<opt_ra_i_r<Z2,Column_types::LIST> >,
							Matrix<opt_ra_i_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z2,Column_types::LIST> >,
							Matrix<opt_ra_r<Z2,Column_types::SET> >,
							Matrix<opt_ra_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i<Z2,Column_types::LIST> >,
							Matrix<opt_ra_i<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z2,Column_types::LIST> >,
							Matrix<opt_ra<Z2,Column_types::SET> >,
							Matrix<opt_ra<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::LIST> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::VECTOR> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_cc_ra<Z2,Column_types::LIST> >,
							Matrix<opt_cc_ra<Z2,Column_types::SET> >,
							Matrix<opt_cc_ra<Z2,Column_types::VECTOR> >,
							Matrix<opt_cc_ra<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_cc_ra<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_cc_ra<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_z2_matrix_types_with_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_base_row_access_columns_options, Matrix, list_of_z2_matrix_types_with_row_access) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	std::vector<std::vector<unsigned int> > rows;
	if constexpr (Matrix::Option_list::has_column_compression){
		//if the union find structure changes, the column_index values of de cells could also change. Change the test with all possibilities?
		rows.push_back({9,10});
		rows.push_back({8,9});
		rows.push_back({8});
		rows.push_back({7});
		rows.push_back({7,9,10});
		rows.push_back({7});
	} else {
		rows.push_back({3,6,9,10,11});
		rows.push_back({3,5,8,9,11});
		rows.push_back({5,8});
		rows.push_back({7});
		rows.push_back({3,6,7,9,10,11});
		rows.push_back({7});
	}

	//rows are unordered
	std::vector<std::set<unsigned int> > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		for (auto& cell : mb.get_row(i)){
			ordered_rows[i].insert(cell.get_column_index());
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		auto& row = ordered_rows[i];
		BOOST_CHECK_EQUAL(row.size(), rows[i].size());
		auto itRow = row.begin();
		auto itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(*itRow, *itChain);
			++itRow; ++itChain;
		}
	}
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Removable_columns_options, Matrix, list_of_matrix_types_with_removable_columns) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 12);
	unsigned int i = 0;
	for (auto& b : ordered_boundaries){
		const auto& col = mb.get_column(i++);
		test_equal_columns(b, col);
	}

	auto it = ordered_boundaries.begin();
	mb.erase_column(10);
	std::advance(it, 10);
	it = ordered_boundaries.erase(it);
	mb.erase_column(7);
	std::advance(it, -3);
	ordered_boundaries.erase(it);
	mb.erase_column(0);
	ordered_boundaries.erase(ordered_boundaries.begin());
	mb.erase_row(5);

	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 9);
	i = 0;
	for (auto& b : ordered_boundaries){
		while (i == 10 || i == 7 || i == 0) ++i;
		const auto& col = mb.get_column(i++);
		test_equal_columns(b, col);
	}

	it = ordered_boundaries.begin();
	mb.erase_column(8);
	std::advance(it, 6);
	it = ordered_boundaries.erase(it);
	mb.erase_column(5);
	std::advance(it, -2);
	ordered_boundaries.erase(it);
	mb.erase_row(2);

	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 7);
	i = 0;
	for (auto& b : ordered_boundaries){
		while (i == 10 || i == 7 || i == 0 || i == 8 || i == 5) ++i;
		const auto& col = mb.get_column(i++);
		test_equal_columns(b, col);
	}
}

//TODO: test removal with row access

