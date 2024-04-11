/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_MATRIX_TESTS_H
#define PM_MATRIX_TESTS_H

#include <stdexcept>
#include <utility>	//std::swap, std::move & std::exchange
// #include <type_traits>
#include <vector>
#include <set>

#include <boost/test/unit_test.hpp>
#include "pm_test_utilities.h"

template<class Column>
std::vector<witness_content<Column> > build_general_matrix(){
	std::vector<witness_content<Column> > columns;

	if constexpr (is_z2<Column>()){
		columns.emplace_back();
		columns.emplace_back();
		columns.emplace_back();
		columns.push_back({0,1,4});
		columns.emplace_back();
		columns.push_back({1,2});
		columns.push_back({0,4});
		columns.push_back({3,4,5});
		columns.push_back({1,2});
		columns.push_back({0,1,4});
		columns.push_back({0,4});
		columns.push_back({0,1,4});
	} else {
		columns.emplace_back();
		columns.emplace_back();
		columns.emplace_back();
		columns.push_back({{0,1},{1,4},{4,1}});
		columns.emplace_back();
		columns.push_back({{1,1},{2,4}});
		columns.push_back({{0,1},{4,4}});
		columns.push_back({{3,1},{4,1},{5,4}});
		columns.push_back({{1,1},{2,4}});
		columns.push_back({{0,1},{1,4},{4,1}});
		columns.push_back({{0,1},{4,4}});
		columns.push_back({{0,1},{1,4},{4,1}});
	}

	return columns;
}

template<class Column>
std::vector<witness_content<Column> > build_simple_boundary_matrix(){
	std::vector<witness_content<Column> > boundaries;

	if constexpr (is_z2<Column>()){
		boundaries.emplace_back();
		boundaries.emplace_back();
		boundaries.emplace_back();
		boundaries.push_back({0,1});
		boundaries.push_back({1,2});
		boundaries.push_back({0,2});
		boundaries.push_back({3,4,5});
	} else {
		// using F = typename Column::Field_element_type;
		// using cont = std::vector<std::pair<unsigned int,F> >;
		// boundaries.emplace_back();
		// boundaries.emplace_back();
		// boundaries.emplace_back();
		// boundaries.push_back(cont{{0,F(1)},{1,F(4)}});
		// boundaries.push_back(cont{{1,F(1)},{2,F(4)}});
		// boundaries.push_back(cont{{0,F(1)},{2,F(4)}});
		// boundaries.push_back(cont{{3,F(1)},{4,F(1)},{5,F(4)}});
		boundaries.emplace_back();
		boundaries.emplace_back();
		boundaries.emplace_back();
		boundaries.push_back({{0,1},{1,4}});
		boundaries.push_back({{1,1},{2,4}});
		boundaries.push_back({{0,1},{2,4}});
		boundaries.push_back({{3,1},{4,1},{5,4}});
	}

	return boundaries;
}

template<class Column>
std::vector<witness_content<Column> > build_simple_chain_matrix(){
	std::vector<witness_content<Column> > columns;

	if constexpr (is_z2<Column>()){
		columns.push_back({0});
		columns.push_back({0,1});
		columns.push_back({0,2});
		columns.push_back({3});
		columns.push_back({3,4});
		columns.push_back({3,4,5});
		columns.push_back({6});
	} else {
		columns.push_back({{0,1}});
		columns.push_back({{0,1},{1,4}});
		columns.push_back({{0,1},{2,4}});
		columns.push_back({{3,1}});
		columns.push_back({{3,1},{4,1}});
		columns.push_back({{3,1},{4,1},{5,4}});
		columns.push_back({{6,1}});
	}

	return columns;
}

template<class Column>
std::vector<witness_content<Column> > build_simple_reduced_row_matrix(){
	std::vector<witness_content<Column> > rows;

	if constexpr (is_z2<Column>()){
		rows.push_back({3});
		rows.push_back({3,4});
		rows.push_back({4});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({6});
	} else {
		rows.push_back({{3,1}});
		rows.push_back({{3,4},{4,1}});
		rows.push_back({{4,4}});
		rows.push_back({{6,1}});
		rows.push_back({{6,1}});
		rows.push_back({{6,4}});
	}

	return rows;
}

template<class Column>
std::vector<witness_content<Column> > build_simple_row_matrix(){
	std::vector<witness_content<Column> > rows;

	if constexpr (is_z2<Column>()){
		rows.push_back({3,5});
		rows.push_back({3,4});
		rows.push_back({4,5});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({6});
	} else {
		rows.push_back({{3,1},{5,1}});
		rows.push_back({{3,4},{4,1}});
		rows.push_back({{4,4},{5,4}});
		rows.push_back({{6,1}});
		rows.push_back({{6,1}});
		rows.push_back({{6,4}});
	}

	return rows;
}

template<class Column>
std::vector<witness_content<Column> > build_simple_chain_row_matrix(){
	std::vector<witness_content<Column> > rows;

	if constexpr (is_z2<Column>()){
		rows.push_back({0,1,2});
		rows.push_back({1});
		rows.push_back({2});
		rows.push_back({3,4,5});
		rows.push_back({4,5});
		rows.push_back({5});
		rows.push_back({6});
	} else {
		rows.push_back({{0,1},{1,1},{2,1}});
		rows.push_back({{1,4}});
		rows.push_back({{2,4}});
		rows.push_back({{3,1},{4,1},{5,1}});
		rows.push_back({{4,1},{5,1}});
		rows.push_back({{5,4}});
		rows.push_back({{6,1}});
	}

	return rows;
}

template<class Column>
std::vector<witness_content<Column> > build_longer_boundary_matrix(){
	std::vector<witness_content<Column> > boundaries;

	if constexpr (is_z2<Column>()){
		boundaries.emplace_back();
		boundaries.emplace_back();
		boundaries.emplace_back();
		boundaries.push_back({0,1});
		boundaries.push_back({1,2});
		boundaries.push_back({0,2});
		boundaries.push_back({3,4,5});
		boundaries.emplace_back();
		boundaries.push_back({1,7});
	} else {
		boundaries.emplace_back();
		boundaries.emplace_back();
		boundaries.emplace_back();
		boundaries.push_back({{0,1},{1,4}});
		boundaries.push_back({{1,1},{2,4}});
		boundaries.push_back({{0,1},{2,4}});
		boundaries.push_back({{3,1},{4,1},{5,4}});
		boundaries.emplace_back();
		boundaries.push_back({{1,1},{7,4}});
	}

	return boundaries;
}

template<class Column>
std::vector<witness_content<Column> > build_longer_chain_matrix(){
	std::vector<witness_content<Column> > columns;

	if constexpr (is_z2<Column>()){
		columns.push_back({0});
		columns.push_back({0,1});
		columns.push_back({0,2});
		columns.push_back({3});
		columns.push_back({3,4});
		columns.push_back({3,4,5});
		columns.push_back({6});
		columns.push_back({0,7});
		columns.push_back({3,8});
	} else {
		columns.push_back({{0,1}});
		columns.push_back({{0,1},{1,4}});
		columns.push_back({{0,1},{2,4}});
		columns.push_back({{3,1}});
		columns.push_back({{3,1},{4,1}});
		columns.push_back({{3,1},{4,1},{5,4}});
		columns.push_back({{6,1}});
		columns.push_back({{0,1},{7,4}});
		columns.push_back({{3,1},{8,1}});
	}

	return columns;
}

template<class Column>
std::vector<witness_content<Column> > build_longer_reduced_row_matrix(){
	std::vector<witness_content<Column> > rows;

	if constexpr (is_z2<Column>()){
		rows.push_back({3});
		rows.push_back({3,4,8});
		rows.push_back({4});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({});
		rows.push_back({8});
	} else {
		rows.push_back({{3,1}});
		rows.push_back({{3,4},{4,1},{8,1}});
		rows.push_back({{4,4}});
		rows.push_back({{6,1}});
		rows.push_back({{6,1}});
		rows.push_back({{6,4}});
		rows.push_back({});
		rows.push_back({{8,4}});
	}

	return rows;
}

template<class Column>
std::vector<witness_content<Column> > build_longer_chain_row_matrix(){
	std::vector<witness_content<Column> > rows;

	if constexpr (is_z2<Column>()){
		rows.push_back({0,1,2,7});
		rows.push_back({1});
		rows.push_back({2});
		rows.push_back({3,4,5,8});
		rows.push_back({4,5});
		rows.push_back({5});
		rows.push_back({6});
		rows.push_back({7});
		rows.push_back({8});
	} else {
		rows.push_back({{0,1},{1,1},{2,1},{7,1}});
		rows.push_back({{1,4}});
		rows.push_back({{2,4}});
		rows.push_back({{3,1},{4,1},{5,1},{8,1}});
		rows.push_back({{4,1},{5,1}});
		rows.push_back({{5,4}});
		rows.push_back({{6,1}});
		rows.push_back({{7,4}});
		rows.push_back({{8,1}});
	}

	return rows;
}

// //*** constructors
// 	Matrix();
// 	template<class Container_type = boundary_type>
// 	Matrix(const std::vector<Container_type>& boundaries);	//simplex indices have to start at 0 and be consecutifs
// 	Matrix(int numberOfColumns);
// 	Matrix(const Matrix &matrixToCopy);
// 	Matrix(Matrix&& other) noexcept;
// 	Matrix& operator=(Matrix other);
// 	friend void swap(Matrix& matrix1, Matrix& matrix2){
// 		swap(matrix1.matrix_, matrix2.matrix_);
// 	}
// 	//******** chain only
// 	template<typename BirthComparatorFunction, typename DeathComparatorFunction>
// 	Matrix(BirthComparatorFunction&& birthComparator, DeathComparatorFunction&& deathComparator = _no_G_death_comparator);
// 	template<typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_type = boundary_type>
// 	Matrix(const std::vector<Boundary_type>& orderedBoundaries,BirthComparatorFunction&& birthComparator, DeathComparatorFunction&& deathComparator = _no_G_death_comparator);
// 	template<typename BirthComparatorFunction, typename DeathComparatorFunction>
// 	Matrix(unsigned int numberOfColumns,BirthComparatorFunction&& birthComparator, DeathComparatorFunction&& deathComparator = _no_G_death_comparator);
// 	//*******************

template<class Matrix>
void test_constructors(){
	std::vector<witness_content<typename Matrix::Column_type> > empty;
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();

	//default constructor
	Matrix m;
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 0);
	test_content_equality(empty, m);

	//constructor from given boundary matrix
	Matrix mb(columns, 5);
	if constexpr (is_RU<Matrix>()){
		columns[5].clear();
	} else if constexpr (is_Chain<Matrix>()){
		columns = build_simple_chain_matrix<typename Matrix::Column_type>();
	}

	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 7);
	test_content_equality(columns, mb);

	//constructor reserving column space
	Matrix mr(5);
	BOOST_CHECK_EQUAL(mr.get_number_of_columns(), 0);
	test_content_equality(empty, mr);

	//copy constructor
	Matrix mc1(mb);
	Matrix mc2 = mb;
	BOOST_CHECK_EQUAL(mc1.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mc2.get_number_of_columns(), 7);
	test_content_equality(columns, mc1);
	test_content_equality(columns, mc2);

	//move constructor
	Matrix mm(std::move(mb));
	BOOST_CHECK_EQUAL(mm.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 0);
	test_content_equality(columns, mm);
	test_content_equality(empty, mb);

	//swap
	swap(mm, mb);
	BOOST_CHECK_EQUAL(mm.get_number_of_columns(), 0);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 7);
	test_content_equality(empty, mm);
	test_content_equality(columns, mb);
}

inline bool birth_comp(unsigned int columnIndex1, unsigned int columnIndex2){ return false; };
inline bool death_comp(unsigned int columnIndex1, unsigned int columnIndex2){ return false; };

template<class Matrix>
void test_chain_constructors(){
	auto ordered_boundaries = build_simple_boundary_matrix<typename Matrix::Column_type>();

	//default constructor
	Matrix m(birth_comp, death_comp);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 0);

	//constructor from given @ref boundarymatrix "boundary matrix"
	Matrix mb(ordered_boundaries, birth_comp, death_comp, 5);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mb.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mb.get_column_dimension(6), 2);

	//constructor reserving column space
	Matrix mr(5, birth_comp, death_comp);
	BOOST_CHECK_EQUAL(mr.get_number_of_columns(), 0);

	//copy constructor
	Matrix mc1(mb);
	Matrix mc2 = mb;
	BOOST_CHECK_EQUAL(mc1.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mc1.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mc1.get_column_dimension(6), 2);
	BOOST_CHECK_EQUAL(mc2.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mc2.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mc2.get_column_dimension(6), 2);

	//move constructor
	Matrix mm(std::move(mb));
	BOOST_CHECK_EQUAL(mm.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mm.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mm.get_column_dimension(6), 2);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 0);

	//swap
	swap(mm, mb);
	BOOST_CHECK_EQUAL(mm.get_number_of_columns(), 0);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mb.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mb.get_column_dimension(6), 2);
}

//*** modifiers & access

// //all options

// 	//base
// 	//base comp
// 	template<class Container_type>
// 	void insert_column(const Container_type& column);
// 	//*******************
// 	//base
// 	//base comp
// 	//boundary
// 	//ru
// 	//chain
// 	//id to pos
// 	//pos to id
// 	unsigned int get_number_of_columns() const;
// 	//base
// 	//base comp
// 	//boundary
// 	//ru: inR = true forced
// 	//chain
// 	//id to pos
// 	//pos to id
// 	bool is_zero_cell(index columnIndex, index rowIndex) const;
// 	//base
// 	//base comp
// 	//boundary
// 	//ru: inR = true forced
// 	//chain: just for sanity checks as a valid @ref chainmatrix "chain matrix" never has an empty column.
// 	//id to pos
// 	//pos to id
// 	bool is_zero_column(index columnIndex);
// 	//*******************

//for base and base comp
template<class Matrix>
void test_general_insertion(){
	auto columns = build_general_matrix<typename Matrix::Column_type>();
	auto col2 = columns.back();
	columns.pop_back();
	auto col1 = columns.back();
	columns.pop_back();

	Matrix m(columns, 5);
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
	m.insert_column(col1);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 11);
	m.insert_column(col2);
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

	// //base: same as insert_column
	// //base comp: same as insert_column
	// //boundary: does not update barcode as it needs reduction
	// //ru
	// //chain: new simplex = new ID even if the same simplex was already inserted and then removed, ie., an ID cannot come back.
	// //id to pos
	// //pos to id
	// template<class Boundary_type = boundary_type>
	// insertion_return_type insert_boundary(const Boundary_type& boundary);
	// //*******************
	// //boundary
	// //ru
	// //chain
	// //id to pos
	// //pos to id
	// dimension_type get_column_dimension(index columnIndex) const;
	// //boundary
	// //ru: only in R
	// //chain
	// //id to pos
	// //pos to id
	// index get_pivot(index columnIndex);
	// //*******************
	// //ru: assumes that pivot exists
	// //chain: assumes that pivot exists
	// //id to pos
	// //pos to id
	// index get_column_with_pivot(index simplexIndex) const;

//for boundary and ru
template<class Matrix>
void test_boundary_insertion(){
	auto orderedBoundaries = build_simple_boundary_matrix<typename Matrix::Column_type>();
	auto boundary2 = orderedBoundaries.back();
	orderedBoundaries.pop_back();
	auto boundary1 = orderedBoundaries.back();
	orderedBoundaries.pop_back();

	Matrix m(orderedBoundaries, 5);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(!m.is_zero_column(4));

	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 5);
	m.insert_boundary(boundary1);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 6);
	m.insert_boundary(boundary2);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 7);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(!m.is_zero_column(4));
	if constexpr (is_RU<Matrix>()){
		BOOST_CHECK(m.is_zero_column(5));	//was reduced
	} else {
		BOOST_CHECK(!m.is_zero_column(5));	//not reduced
	}
	BOOST_CHECK(!m.is_zero_column(6));

	BOOST_CHECK_EQUAL(m.get_column_dimension(0), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(1), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(3), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(4), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(5), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(6), 2);

	BOOST_CHECK_EQUAL(m.get_pivot(0), static_cast<typename Matrix::id_index>(-1));
	BOOST_CHECK_EQUAL(m.get_pivot(1), static_cast<typename Matrix::id_index>(-1));
	BOOST_CHECK_EQUAL(m.get_pivot(2), static_cast<typename Matrix::id_index>(-1));
	BOOST_CHECK_EQUAL(m.get_pivot(3), 1);
	BOOST_CHECK_EQUAL(m.get_pivot(4), 2);
	if constexpr (is_RU<Matrix>()){
		BOOST_CHECK_EQUAL(m.get_pivot(5), static_cast<typename Matrix::id_index>(-1));	//was reduced
	} else {
		BOOST_CHECK_EQUAL(m.get_pivot(5), 2);	//not reduced
	}
	BOOST_CHECK_EQUAL(m.get_pivot(6), 5);

	if constexpr (is_RU<Matrix>()){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 6);
	}
}

	// //chain: new simplex = new ID even if the same simplex was already inserted and then removed, ie., an ID cannot come back.
	// //id to pos
	// template<class Boundary_type>
	// insertion_return_type insert_boundary(index simplexIndex, const Boundary_type& boundary);
	// //*******************

//for chain
template<class Matrix>
void test_chain_boundary_insertion(Matrix& m1, Matrix& m2){
	auto test = [](Matrix& m){
		BOOST_CHECK(m.is_zero_cell(1,2));
		BOOST_CHECK(!m.is_zero_cell(1,1));
		BOOST_CHECK(!m.is_zero_cell(3,3));
		BOOST_CHECK(!m.is_zero_column(0));
		BOOST_CHECK(!m.is_zero_column(1));
		BOOST_CHECK(!m.is_zero_column(2));
		BOOST_CHECK(!m.is_zero_column(3));
		BOOST_CHECK(!m.is_zero_column(4));
		BOOST_CHECK(!m.is_zero_column(5));
		BOOST_CHECK(!m.is_zero_column(6));

		BOOST_CHECK_EQUAL(m.get_column_dimension(0), 0);
		BOOST_CHECK_EQUAL(m.get_column_dimension(1), 0);
		BOOST_CHECK_EQUAL(m.get_column_dimension(2), 0);
		BOOST_CHECK_EQUAL(m.get_column_dimension(3), 1);
		BOOST_CHECK_EQUAL(m.get_column_dimension(4), 1);
		BOOST_CHECK_EQUAL(m.get_column_dimension(5), 1);
		BOOST_CHECK_EQUAL(m.get_column_dimension(6), 2);

		BOOST_CHECK_EQUAL(m.get_pivot(0), 0);
		BOOST_CHECK_EQUAL(m.get_pivot(1), 1);
		BOOST_CHECK_EQUAL(m.get_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_pivot(3), 3);
		BOOST_CHECK_EQUAL(m.get_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_pivot(5), 5);
		BOOST_CHECK_EQUAL(m.get_pivot(6), 6);

		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
	};

	auto orderedBoundaries = build_simple_boundary_matrix<typename Matrix::Column_type>();

	for (unsigned int i = 0; i < orderedBoundaries.size(); ++i){
		if constexpr (is_indexed_by_position<Matrix>()) m1.insert_boundary(orderedBoundaries[i]);
		else m1.insert_boundary(i, orderedBoundaries[i]);
	}
	BOOST_CHECK_EQUAL(m1.get_number_of_columns(), 7);
	test(m1);

	auto boundary2 = orderedBoundaries.back();
	orderedBoundaries.pop_back();
	auto boundary1 = orderedBoundaries.back();
	orderedBoundaries.pop_back();

	BOOST_CHECK(m2.is_zero_cell(1,2));
	BOOST_CHECK(!m2.is_zero_cell(1,1));
	BOOST_CHECK(!m2.is_zero_cell(3,3));
	BOOST_CHECK(!m2.is_zero_column(0));
	BOOST_CHECK(!m2.is_zero_column(1));
	BOOST_CHECK(!m2.is_zero_column(2));
	BOOST_CHECK(!m2.is_zero_column(3));
	BOOST_CHECK(!m2.is_zero_column(4));

	BOOST_CHECK_EQUAL(m2.get_number_of_columns(), 5);
	m2.insert_boundary(boundary1);
	BOOST_CHECK_EQUAL(m2.get_number_of_columns(), 6);
	m2.insert_boundary(boundary2);
	BOOST_CHECK_EQUAL(m2.get_number_of_columns(), 7);
	test(m2);
}

	// //base
	// //base comp: non const because of path compression in union-find
	// //boundary
	// //ru: inR = true forced
	// //chain
	// //id to pos
	// //pos to id
	// returned_column_type& get_column(index columnIndex);
	// //*******************
	// //base
	// //boundary
	// //ru: inR = true forced
	// //chain
	// //id to pos
	// //pos to id
	// const Column_type& get_column(index columnIndex) const;
	// //*******************

template<class Matrix>
void test_base_access(){
	auto columns = build_general_matrix<typename Matrix::Column_type>();
	auto col2 = columns.back();
	columns.pop_back();
	auto col1 = columns.back();
	columns.pop_back();

	Matrix m(columns, 5);

	test_content_equality(columns, m);

	m.insert_column(col1);
	columns.push_back(col1);
	m.insert_column(col2);
	columns.push_back(col2);

	test_content_equality(columns, m);
}

template<class Matrix>
void test_boundary_access(){
	auto orderedBoundaries = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(orderedBoundaries, 5);

	for (unsigned int i = 0; i < orderedBoundaries.size(); ++i){
		if constexpr (is_RU<Matrix>()){
			if (i == 5){
				BOOST_CHECK(m.get_column(i).is_empty());	//reduced
			} else {
				const auto& col = m.get_column(i);	//to force the const version
				test_column_equality<typename Matrix::Column_type>(orderedBoundaries[i], get_column_content_via_iterators(col));
			}
		} else {
			const auto& col = m.get_column(i);	//to force the const version
			test_column_equality<typename Matrix::Column_type>(orderedBoundaries[i], get_column_content_via_iterators(col));
		}
	}
}

template<class Matrix>
void test_chain_access(Matrix& m){
	auto columns = build_simple_chain_matrix<typename Matrix::Column_type>();
	test_content_equality(columns, m);
}

	// //base
	// //boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	// //ru: inR = true forced, avoid calling with specialized options or make it such that it makes sense for persistence
	// //id to pos
	// void zero_cell(index columnIndex, index rowIndex);
	// //base
	// //boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	// //ru: inR = true forced, avoid calling with specialized options or make it such that it makes sense for persistence
	// //id to pos
	// void zero_column(index columnIndex);
	// //*******************

template<class Matrix>
void test_zeroing(){
	auto orderedBoundaries = build_simple_boundary_matrix<typename Matrix::Column_type>();

	Matrix m(orderedBoundaries, 5);

	BOOST_CHECK(!m.is_zero_cell(3,1));
	BOOST_CHECK(!m.is_zero_column(3));
	m.zero_cell(3, 1);
	BOOST_CHECK(m.is_zero_cell(3,1));
	BOOST_CHECK(!m.is_zero_column(3));
	m.zero_column(3);
	BOOST_CHECK(m.is_zero_cell(3,1));
	BOOST_CHECK(m.is_zero_column(3));
}

	// //ru
	// Column_type& get_column(index columnIndex, bool inR);
	// //ru
	// const Column_type& get_column(index columnIndex, bool inR) const;
	// //ru
	// void zero_cell(index columnIndex, index rowIndex, bool inR);
	// //ru
	// void zero_column(index columnIndex, bool inR);
	// //ru
	// bool is_zero_cell(index columnIndex, index rowIndex, bool inR) const;
	// //ru
	// bool is_zero_column(index columnIndex, bool inR);
	// //*******************

template<class Matrix>
void test_ru_u_access(){
	auto orderedBoundaries = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(orderedBoundaries, 5);

	std::vector<witness_content<typename Matrix::Column_type> > uColumns(7);
	if constexpr (Matrix::Option_list::is_z2){
		if constexpr (Matrix::Option_list::has_vine_update){
			uColumns[0] = {0};
			uColumns[1] = {1};
			uColumns[2] = {2};
			uColumns[3] = {3,5};
			uColumns[4] = {4,5};
			uColumns[5] = {5};
			uColumns[6] = {6};
		} else {
			uColumns[0] = {0};
			uColumns[1] = {1};
			uColumns[2] = {2};
			uColumns[3] = {3};
			uColumns[4] = {4};
			uColumns[5] = {3,4,5};
			uColumns[6] = {6};
		}
	} else {
		uColumns[0] = {{0,1}};
		uColumns[1] = {{1,1}};
		uColumns[2] = {{2,1}};
		uColumns[3] = {{3,1}};
		uColumns[4] = {{4,1}};
		uColumns[5] = {{3,1},{4,1},{5,4}};
		uColumns[6] = {{6,1}};
	}

	unsigned int i = 0;
	for (auto& c : uColumns){
		const auto& col = m.get_column(i++, false);	//to force the const version
		test_column_equality<typename Matrix::Column_type>(c, get_column_content_via_iterators(col));
	}

	if constexpr (Matrix::Option_list::has_vine_update){
		BOOST_CHECK(!m.is_zero_cell(3, 5, false));
		BOOST_CHECK(!m.is_zero_column(4, false));
		m.zero_cell(3, 5, false);
		BOOST_CHECK(m.is_zero_cell(3, 5, false));
		BOOST_CHECK(!m.is_zero_column(4, false));
		m.zero_column(4, false);
		BOOST_CHECK(m.is_zero_cell(3, 5, false));
		BOOST_CHECK(m.is_zero_column(4, false));
	} else {
		BOOST_CHECK(!m.is_zero_cell(5, 3, false));
		BOOST_CHECK(!m.is_zero_column(4, false));
		m.zero_cell(5, 3, false);
		BOOST_CHECK(m.is_zero_cell(5, 3, false));
		BOOST_CHECK(!m.is_zero_column(4, false));
		m.zero_column(4, false);
		BOOST_CHECK(m.is_zero_cell(5, 3, false));
		BOOST_CHECK(m.is_zero_column(4, false));
	}
}

// //row access
	// 	//base
	// 	//base comp
	// 	//boundary
	// 	//ru: inR = true forced
	// 	//chain
	// 	//id to pos
	// 	//pos to id
	// 	const Row_type& get_row(index rowIndex) const;
	// 	//*******************
	// 	//base
	// 	//boundary
	// 	//ru: inR = true forced
	// 	//chain
	// 	//id to pos
	// 	//pos to id
	// 	Row_type& get_row(index rowIndex);
	// 	//*******************

template<class Matrix>
void test_base_z2_row_access(){
	auto columns = build_general_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	std::vector<std::vector<unsigned int> > rows;
	if constexpr (Matrix::Option_list::has_column_compression){
		//if the union find structure changes, the column_index values of de cells could also change. Change the test with all possibilities?
		rows.push_back({3,6});
		rows.push_back({3,5});
		rows.push_back({5});
		rows.push_back({7});
		rows.push_back({3,6,7});
		rows.push_back({7});
	} else {
		rows.push_back({3,6,9,10,11});
		rows.push_back({3,5,8,9,11});
		rows.push_back({5,8});
		rows.push_back({7});
		rows.push_back({3,6,7,9,10,11});
		rows.push_back({7});
	}

	unsigned int i = 0;
	for (auto& r : rows){
		test_column_equality<typename Matrix::Column_type>(r, get_ordered_row(m, i++));
	}
}

template<class Matrix>
void test_base_z5_row_access(){
	auto columns = build_general_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	std::vector<std::vector<std::pair<unsigned int,typename Matrix::element_type> > > rows;
	if constexpr (Matrix::Option_list::has_column_compression){
		//if the union find structure changes, the column_index values of de cells could also change. Change the test with all possibilities?
		rows.push_back({{3,1},{6,1}});
		rows.push_back({{3,4},{5,1}});
		rows.push_back({{5,4}});
		rows.push_back({{7,1}});
		rows.push_back({{3,1},{6,4},{7,1}});
		rows.push_back({{7,4}});
	} else {
		rows.push_back({{3,1},{6,1},{9,1},{10,1},{11,1}});
		rows.push_back({{3,4},{5,1},{8,1},{9,4},{11,4}});
		rows.push_back({{5,4},{8,4}});
		rows.push_back({{7,1}});
		rows.push_back({{3,1},{6,4},{7,1},{9,1},{10,4},{11,1}});
		rows.push_back({{7,4}});
	}

	unsigned int i = 0;
	for (auto& r : rows){
		test_column_equality<typename Matrix::Column_type>(r, get_ordered_row(m, i++));
	}
}

template<class Matrix>
void test_non_base_row_access(Matrix& m){
	std::vector<witness_content<typename Matrix::Column_type> > rows;
	if constexpr (Matrix::Option_list::is_of_boundary_type){
		if constexpr (is_RU<Matrix>()){
			rows = build_simple_reduced_row_matrix<typename Matrix::Column_type>();
		} else {
			rows = build_simple_row_matrix<typename Matrix::Column_type>();
		}
	} else {
		rows = build_simple_chain_row_matrix<typename Matrix::Column_type>();
	}
	

	unsigned int i = 0;
	for (auto& r : rows){
		test_column_equality<typename Matrix::Column_type>(r, get_ordered_row(m, i++));
	}
}

	// 	//ru
	// 	Row_type& get_row(index rowIndex, bool inR);
	// 	//ru
	// 	const Row_type& get_row(index rowIndex, bool inR) const;
	// 	//*******************

template<class Matrix>
void test_ru_u_row_access(){
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	std::vector<witness_content<typename Matrix::Column_type> > rows;
	if constexpr (Matrix::Option_list::is_z2){
		if constexpr (Matrix::Option_list::has_vine_update){
			rows.push_back({0});
			rows.push_back({1});
			rows.push_back({2});
			rows.push_back({3});
			rows.push_back({4});
			rows.push_back({3,4,5});
			rows.push_back({6});
		} else {
			rows.push_back({0});
			rows.push_back({1});
			rows.push_back({2});
			rows.push_back({3,5});
			rows.push_back({4,5});
			rows.push_back({5});
			rows.push_back({6});
		}
	} else {
		rows.push_back({{0,1}});
		rows.push_back({{1,1}});
		rows.push_back({{2,1}});
		rows.push_back({{3,1},{5,1}});
		rows.push_back({{4,1},{5,1}});
		rows.push_back({{5,4}});
		rows.push_back({{6,1}});
	}

	column_content<typename Matrix::Column_type> orderedRows;
	unsigned int i = 0;
	for (auto& r : rows){
		orderedRows.clear();
		for (const auto& cell : m.get_row(i++, false)){
			if constexpr (Matrix::Option_list::is_z2){
				orderedRows.insert(cell.get_column_index());
			} else {
				orderedRows.insert({cell.get_column_index(), cell.get_element()});
			}
		}
		test_column_equality<typename Matrix::Column_type>(r, orderedRows);
	}
}

// //row access + erasable rows

	// 	//base: as base comp but also has to test with swaps
	// 	//base comp: no row access necessary (but does nothing otherwise), assumes the row is empty, just thought as an index cleanup
	// 	//boundary: indirect, assumes the row is empty
	// 	//ru: indirect, assumes the row is empty
	// 	//chain: indirect, assumes the row is empty
	// 	//id to pos
	// 	//pos to id
	// 	void erase_row(index rowIndex);
	// 	//*******************

template<class Matrix>
void test_row_removal(){
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	columns[6].pop_back();	//empties row 5. Not a legit @ref boundarymatrix "boundary matrix" anymore, but for the test, should be fine, except for chain.

	Matrix m(columns, 5);

	m.erase_row(5);

	BOOST_CHECK_THROW(m.get_row(5), std::logic_error);
}

template<class Matrix>
void test_chain_row_removal(Matrix& m){
	m.erase_row(6);	//not empty, so ignored
	BOOST_CHECK_NO_THROW(m.get_row(6));

	if constexpr (Matrix::Option_list::has_map_column_container || !Matrix::Option_list::has_vine_update){
		m.remove_last();	//calls erase_row(6)
		BOOST_CHECK_THROW(m.get_row(6), std::logic_error);
	}
}

// //erasable columns

	// 	//base
	// 	template<class Container_type>
	// 	void insert_column(const Container_type& column, int columnIndex);
	// 	//*******************
	// 	//base: (as direct link between column and row is not guaranteed)
	// 	void remove_column(index columnIndex);
	// 	//*******************

template<class Matrix>
void test_column_removal(){
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	test_content_equality(columns, m);

	m.remove_column(2);
	m.remove_column(4);
	BOOST_CHECK_THROW(m.get_column(2), std::logic_error);
	BOOST_CHECK_THROW(m.get_column(4), std::logic_error);

	columns.erase(columns.begin() + 4);
	columns.erase(columns.begin() + 2);
	unsigned int i = 0;
	for (auto& b : columns){
		if (i == 2 || i == 4) ++i;
		test_column_equality<typename Matrix::Column_type>(b, get_column_content_via_iterators(m.get_column(i++)));
	}

	if constexpr (!Matrix::Option_list::has_row_access){
		m.insert_column(witness_content<typename Matrix::Column_type>{}, 2);
		BOOST_CHECK_NO_THROW(m.get_column(2));
		BOOST_CHECK_THROW(m.get_column(4), std::logic_error);
	}
}

	// //boundary: update barcode if already computed, does not verify if it really was maximal
	// //ru
	// //chain
	// //id to pos
	// //pos to id
	// void remove_maximal_face(index columnIndex);
	// //*******************

template<class Matrix>
void test_boundary_maximal_simplex_removal(){
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	test_content_equality(columns, m);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(m.get_current_barcode().back().death, 6);	//pairing always true for boundary for now (only thing differenciating it from base)

	m.remove_last();

	columns.pop_back();
	columns[5].clear();	//was reduced

	test_content_equality(columns, m);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 6);
	BOOST_CHECK_EQUAL(m.get_current_barcode().back().death, static_cast<typename Matrix::pos_index>(-1));
}

template<class Matrix>
void test_ru_maximal_simplex_removal(){
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	columns[5].clear();

	test_content_equality(columns, m);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 7);
	if constexpr (Matrix::Option_list::has_column_pairings){
		BOOST_CHECK_EQUAL(m.get_current_barcode().back().death, 6);
	}

	if constexpr (Matrix::Option_list::has_vine_update){
		m.remove_maximal_face(6);
	} else {
		m.remove_last();
	}

	columns.pop_back();
	test_content_equality(columns, m);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 6);
	if constexpr (Matrix::Option_list::has_column_pairings){
		BOOST_CHECK_EQUAL(m.get_current_barcode().back().death, static_cast<typename Matrix::pos_index>(-1));
	}
}

template<class Matrix>
void test_chain_maximal_simplex_removal(Matrix& m){
	auto columns = build_simple_chain_matrix<typename Matrix::Column_type>();

	test_content_equality(columns, m);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 7);
	if constexpr (Matrix::Option_list::has_column_pairings){
		BOOST_CHECK_EQUAL(m.get_current_barcode().back().death, 6);
	}

	if constexpr (Matrix::Option_list::has_vine_update && Matrix::Option_list::has_map_column_container && Matrix::Option_list::has_column_pairings){
		m.remove_maximal_face(6);
	} else {
		m.remove_last();
	}

	columns.pop_back();
	test_content_equality(columns, m);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 6);
	if constexpr (Matrix::Option_list::has_column_pairings){
		BOOST_CHECK_EQUAL(m.get_current_barcode().back().death, static_cast<typename Matrix::pos_index>(-1));
	}
}

// //max dimension

	// //boundary: indirect
	// //ru
	// //chain: indirect
	// //id to pos
	// //pos to id
	// dimension_type get_max_dimension() const;
	// //*******************

template<class Matrix>
void test_maximal_dimension(Matrix& m){
	BOOST_CHECK_EQUAL(m.get_max_dimension(), 2);

	if constexpr (Matrix::Option_list::is_z2){
		m.insert_boundary({0,1,2,3,4,5});
		m.insert_boundary({0,1,2,3});
	} else {
		m.insert_boundary({{0,1},{1,4},{2,1},{3,4},{4,1},{5,4}});
		m.insert_boundary({{0,1},{1,4},{2,1},{3,4}});
	}

	BOOST_CHECK_EQUAL(m.get_max_dimension(), 5);

	if constexpr (Matrix::Option_list::has_vine_update && (Matrix::Option_list::is_of_boundary_type || (Matrix::Option_list::has_map_column_container && Matrix::Option_list::has_column_pairings))){
		m.remove_maximal_face(7);
		BOOST_CHECK_EQUAL(m.get_max_dimension(), 3);
	}
}

//*** operators

	// //base
	// //base comp
	// //boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	// //ru: avoid calling with specialized options or make it such that it makes sense for persistence
	// //chain: avoid calling with specialized options or make it such that it makes sense for persistence
	// //id to pos
	// //pos to id
	// template<typename Index_type>
	// std::enable_if_t<std::is_integral_v<Index_type> > add_to(Index_type sourceColumnIndex, Index_type targetColumnIndex);
	// //base: necessary because of vector columns ?
	// //base comp: necessary because of vector columns ?
	// //boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	// //ru: avoid calling with specialized options or make it such that it makes sense for persistence
	// //chain: avoid calling with specialized options or make it such that it makes sense for persistence
	// //id to pos:indirect
	// //pos to id
	// void add_to(Column_type& sourceColumn, index targetColumnIndex);
	// //base
	// //base comp
	// //boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	// //ru: avoid calling with specialized options or make it such that it makes sense for persistence
	// //chain: avoid calling with specialized options or make it such that it makes sense for persistence
	// //id to pos: indirect
	// //pos to id
	// void add_to(Column_type& sourceColumn, const Field_type& coefficient, index targetColumnIndex);
	// //base
	// //base comp
	// //boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	// //ru: avoid calling with specialized options or make it such that it makes sense for persistence
	// //chain: avoid calling with specialized options or make it such that it makes sense for persistence
	// //id to pos: indirect
	// //pos to id
	// void add_to(const Field_type& coefficient, Column_type& sourceColumn, index targetColumnIndex);
	//*******************

template<class Matrix>
void test_base_operation(){
	auto columns = build_general_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	test_content_equality(columns, m);

	m.add_to(6, 7);
	if constexpr (Matrix::Option_list::is_z2){
		columns[7] = {0,3,5};
	} else {
		columns[7] = {{0,1},{3,1},{5,4}};
	}
	test_content_equality(columns, m);

	if constexpr (!Matrix::isNonBasic){
		m.add_to(m.get_column(3), 10);
		if constexpr (Matrix::Option_list::is_z2){
			columns[10] = {1};
		} else {
			columns[10] = {{0,2},{1,4}};
		}
		test_content_equality(columns, m);
	}

	m.multiply_target_and_add_to(3, 3, 5);
	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {0,2,4};
	} else {
		columns[5] = {{0,1},{1,2},{2,2},{4,1}};
	}
	test_content_equality(columns, m);

	m.multiply_source_and_add_to(3, 3, 6);
	if constexpr (Matrix::Option_list::is_z2){
		columns[6] = {1};
	} else {
		columns[6] = {{0,4},{1,2},{4,2}};
	}
	test_content_equality(columns, m);
}

template<class Matrix>
void test_base_col_comp_operation(){
	auto columns = build_general_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	test_content_equality(columns, m);

	m.add_to(6, 7);
	if constexpr (Matrix::Option_list::is_z2){
		columns[7] = {0,3,5};
	} else {
		columns[7] = {{0,1},{3,1},{5,4}};
	}
	test_content_equality(columns, m);

	m.add_to(m.get_column(3), 10);
	if constexpr (Matrix::Option_list::is_z2){
		columns[6] = {1};
		columns[10] = {1};
	} else {
		columns[6] = {{0,2},{1,4}};
		columns[10] = {{0,2},{1,4}};
	}
	test_content_equality(columns, m);

	m.multiply_target_and_add_to(m.get_column(3), 3, 5);
	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {0,2,4};
		columns[8] = {0,2,4};
	} else {
		columns[5] = {{0,1},{1,2},{2,2},{4,1}};
		columns[8] = {{0,1},{1,2},{2,2},{4,1}};
	}
	test_content_equality(columns, m);

	m.multiply_source_and_add_to(3, m.get_column(3), 6);
	if constexpr (Matrix::Option_list::is_z2){
		columns[6] = {0,4};
		columns[10] = {0,4};
	} else {
		columns[6] = {{1,1},{4,3}};
		columns[10] = {{1,1},{4,3}};
	}
	test_content_equality(columns, m);
}

template<class Matrix>
void test_ru_operation(){
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	columns[5].clear();

	std::vector<witness_content<typename Matrix::Column_type> > uColumns(7);
	if constexpr (Matrix::Option_list::is_z2){
		if constexpr (Matrix::Option_list::has_vine_update){
			uColumns[0] = {0};
			uColumns[1] = {1};
			uColumns[2] = {2};
			uColumns[3] = {3,5};
			uColumns[4] = {4,5};
			uColumns[5] = {5};
			uColumns[6] = {6};
		} else {
			uColumns[0] = {0};
			uColumns[1] = {1};
			uColumns[2] = {2};
			uColumns[3] = {3};
			uColumns[4] = {4};
			uColumns[5] = {3,4,5};
			uColumns[6] = {6};
		}
	} else {
		uColumns[0] = {{0,1}};
		uColumns[1] = {{1,1}};
		uColumns[2] = {{2,1}};
		uColumns[3] = {{3,1}};
		uColumns[4] = {{4,1}};
		uColumns[5] = {{3,1},{4,1},{5,4}};
		uColumns[6] = {{6,1}};
	}

	test_content_equality(columns, m);
	unsigned int i = 0;
	if constexpr (is_indexed_by_position<Matrix>()){
		for (auto& b : uColumns){
			test_column_equality<typename Matrix::Column_type>(b, get_column_content_via_iterators(m.get_column(i++, false)));
		}
	}

	m.add_to(3, 5);
	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {0,1};
		if constexpr (Matrix::Option_list::has_vine_update)
			uColumns[3] = {3};
		else
			uColumns[5] = {4,5};
	} else {
		columns[5] = {{0,1},{1,4}};
		uColumns[5] = {{3,2},{4,1},{5,4}};
	}
	test_content_equality(columns, m);
	if constexpr (is_indexed_by_position<Matrix>()){
		i = 0;
		for (auto& b : uColumns){
			test_column_equality<typename Matrix::Column_type>(b, get_column_content_via_iterators(m.get_column(i++, false)));
		}
	}

	m.add_to(4, 5);
	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {0,2};
		if constexpr (Matrix::Option_list::has_vine_update)
			uColumns[4] = {4};
		else
			uColumns[5] = {5};
	} else {
		columns[5] = {{0,1},{2,4}};
		uColumns[5] = {{3,2},{4,2},{5,4}};
	}
	test_content_equality(columns, m);
	if constexpr (is_indexed_by_position<Matrix>()){
		i = 0;
		for (auto& b : uColumns){
			test_column_equality<typename Matrix::Column_type>(b, get_column_content_via_iterators(m.get_column(i++, false)));
		}
	}

	if constexpr (!Matrix::Option_list::has_vine_update){
		m.multiply_target_and_add_to(5, 3, 3);
		if constexpr (Matrix::Option_list::is_z2){
			columns[3] = {1,2};
			uColumns[3] = {3,5};
		} else {
			columns[3] = {{0,4},{1,2},{2,4}};
			uColumns[3] = {{4,2},{5,4}};
		}
		test_content_equality(columns, m);
		if constexpr (is_indexed_by_position<Matrix>()){
			i = 0;
			for (auto& b : uColumns){
				test_column_equality<typename Matrix::Column_type>(b, get_column_content_via_iterators(m.get_column(i++, false)));
			}
		}

		m.multiply_source_and_add_to(4, 5, 4);
		if constexpr (Matrix::Option_list::is_z2){
			columns[4] = {1,2};
			uColumns[4] = {4};
		} else {
			columns[4] = {{0,4},{1,1}};
			uColumns[4] = {{3,3},{4,4},{5,1}};
		}
		test_content_equality(columns, m);
		if constexpr (is_indexed_by_position<Matrix>()){
			i = 0;
			for (auto& b : uColumns){
				test_column_equality<typename Matrix::Column_type>(b, get_column_content_via_iterators(m.get_column(i++, false)));
			}
		}
	}
}

template<class Matrix>
void test_chain_operation(Matrix& m){
	auto columns = build_simple_chain_matrix<typename Matrix::Column_type>();

	test_content_equality(columns, m);

	m.add_to(3, 5);
	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {4,5};
	} else {
		columns[5] = {{3,2},{4,1},{5,4}};
	}
	test_content_equality(columns, m);

	m.add_to(4, 5);
	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {3,5};
	} else {
		columns[5] = {{3,3},{4,2},{5,4}};
	}
	test_content_equality(columns, m);

	m.multiply_target_and_add_to(5, 3, 3);
	if constexpr (Matrix::Option_list::is_z2){
		columns[3] = {5};
	} else {
		columns[3] = {{3,1},{4,2},{5,4}};
	}
	test_content_equality(columns, m);

	m.multiply_source_and_add_to(4, 5, 4);
	if constexpr (Matrix::Option_list::is_z2){
		columns[4] = {3,4};
	} else {
		columns[4] = {{3,3},{4,4},{5,1}};
	}
	test_content_equality(columns, m);
}

	// //base
	// //base comp
	// template<class Cell_range>
	// std::enable_if_t<!std::is_integral_v<Cell_range> > add_to(const Cell_range& sourceColumn, index targetColumnIndex);
	// //base
	// //base comp
	// template<class Cell_range>
	// void add_to(const Cell_range& sourceColumn, const Field_type& coefficient, index targetColumnIndex);
	// //base
	// //base comp
	// template<class Cell_range>
	// void add_to(const Field_type& coefficient, const Cell_range& sourceColumn, index targetColumnIndex);
	// //*******************

template<class Matrix>
void test_base_cell_range_operation(){
	using Cell = typename Matrix::Cell_type;

	auto columns = build_general_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	std::vector<Cell> range;
	range = {Cell(0),Cell(1),Cell(4)};
	if constexpr (!Matrix::Option_list::is_z2){
		range[0].set_element(1);
		range[1].set_element(4);
		range[2].set_element(1);
	}

	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[10] = {1};
	} else {
		columns[10] = {{0,2},{1,4}};
	}
	m.add_to(range, 10);
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {0,2,4};
	} else {
		columns[5] = {{0,1},{1,2},{2,2},{4,1}};
	}
	m.multiply_target_and_add_to(range, 3, 5);
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[6] = {1};
	} else {
		columns[6] = {{0,4},{1,2},{4,2}};
	}
	m.multiply_source_and_add_to(3, range, 6);
	test_content_equality(columns, m);
}

template<class Matrix>
void test_base_col_comp_cell_range_operation(){
	using Cell = typename Matrix::Cell_type;

	auto columns = build_general_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	std::vector<Cell> range;
	range = {Cell(0),Cell(1),Cell(4)};
	if constexpr (!Matrix::Option_list::is_z2){
		range[0].set_element(1);
		range[1].set_element(4);
		range[2].set_element(1);
	}

	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[6] = {1};
		columns[10] = {1};
	} else {
		columns[6] = {{0,2},{1,4}};
		columns[10] = {{0,2},{1,4}};
	}
	m.add_to(range, 10);
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {0,2,4};
		columns[8] = {0,2,4};
	} else {
		columns[5] = {{0,1},{1,2},{2,2},{4,1}};
		columns[8] = {{0,1},{1,2},{2,2},{4,1}};
	}
	m.multiply_target_and_add_to(range, 3, 5);
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[6] = {0,4};
		columns[10] = {0,4};
	} else {
		columns[6] = {{1,1},{4,3}};
		columns[10] = {{1,1},{4,3}};
	}
	m.multiply_source_and_add_to(3, range, 6);
	test_content_equality(columns, m);
}

	// //base: indirect via Cell_range version
	// //base comp: indirect via Cell_range version
	// //boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	// //ru: avoid calling with specialized options or make it such that it makes sense for persistence
	// //id to pos
	// void add_to(const Column_type& sourceColumn, index targetColumnIndex);
	// //base: indirect via Cell_range version
	// //base comp: indirect via Cell_range version
	// //boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	// //ru: avoid calling with specialized options or make it such that it makes sense for persistence
	// //id to pos
	// void add_to(const Column_type& sourceColumn, const Field_type& coefficient, index targetColumnIndex);
	// //base: indirect via Cell_range version
	// //base comp: indirect via Cell_range version
	// //boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	// //ru: avoid calling with specialized options or make it such that it makes sense for persistence
	// //id to pos
	// void add_to(const Field_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex);

template<class Matrix>
void test_const_operation(){
	using C = typename Matrix::Column_type;
	typename Matrix::Field_operators op(5);
	typename Matrix::Cell_constructor pool;

	auto columns = build_general_matrix<C>();
	Matrix m(columns, 5);

	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[10] = {1};
		m.add_to(C({0,1,4}, nullptr, &pool), 10);	//only works with the const version because of reference
	} else {
		columns[10] = {{0,2},{1,4}};
		m.add_to(C({{0,1},{1,4},{4,1}}, &op, &pool), 10);
	}
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {0,2,4};
		m.multiply_target_and_add_to(C({0,1,4}, nullptr, &pool), 3, 5);
	} else {
		columns[5] = {{0,1},{1,2},{2,2},{4,1}};
		m.multiply_target_and_add_to(C({{0,1},{1,4},{4,1}}, &op, &pool), 3, 5);
	}
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[6] = {1};
		m.multiply_source_and_add_to(3, C({0,1,4}, nullptr, &pool), 6);
	} else {
		columns[6] = {{0,4},{1,2},{4,2}};
		m.multiply_source_and_add_to(3, C({{0,1},{1,4},{4,1}}, &op, &pool), 6);
	}
	test_content_equality(columns, m);
}

template<class Matrix>
void test_base_col_comp_const_operation(){
	using C = typename Matrix::Column_type;
	typename Matrix::Field_operators op(5);
	typename Matrix::Cell_constructor pool;

	auto columns = build_general_matrix<C>();
	Matrix m(columns, 5);

	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[6] = {1};
		columns[10] = {1};
		m.add_to(C({0,1,4}, nullptr, &pool), 10);	//only works with the const version because of reference
	} else {
		columns[6] = {{0,2},{1,4}};
		columns[10] = {{0,2},{1,4}};
		m.add_to(C({{0,1},{1,4},{4,1}}, &op, &pool), 10);
	}
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[5] = {0,2,4};
		columns[8] = {0,2,4};
		m.multiply_target_and_add_to(C({0,1,4}, nullptr, &pool), 3, 5);
	} else {
		columns[5] = {{0,1},{1,2},{2,2},{4,1}};
		columns[8] = {{0,1},{1,2},{2,2},{4,1}};
		m.multiply_target_and_add_to(C({{0,1},{1,4},{4,1}}, &op, &pool), 3, 5);
	}
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::is_z2){
		columns[6] = {0,4};
		columns[10] = {0,4};
		m.multiply_source_and_add_to(3, C({0,1,4}, nullptr, &pool), 6);
	} else {
		columns[6] = {{1,1},{4,3}};
		columns[10] = {{1,1},{4,3}};
		m.multiply_source_and_add_to(3, C({{0,1},{1,4},{4,1}}, &op, &pool), 6);
	}
	test_content_equality(columns, m);
}

// //*** Persistence diagram
// 	//boundary
// 	//ru: indirect through const version
// 	//chain: indirect through const version
// 	//id to pos
// 	//pos to id: indirect
// 	const barcode_type& get_current_barcode();
// 	//*******************
// 	//chain
// 	//ru
// 	//pos to id
// 	const barcode_type& get_current_barcode() const;

template<class Matrix>
void test_barcode(){
	struct BarComp {
		bool operator()(const std::tuple<int,int,int>& c1, const std::tuple<int,int,int>& c2) const
		{
			if (std::get<0>(c1) == std::get<0>(c2))
				return std::get<1>(c1) < std::get<1>(c2);
			return std::get<0>(c1) < std::get<0>(c2);
		}
	};

	auto columns = build_longer_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	const auto& barcode = m.get_current_barcode();

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		columns[5].clear();
	} else {
		columns = build_longer_chain_matrix<typename Matrix::Column_type>();
	}
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::has_row_access){
		std::vector<witness_content<typename Matrix::Column_type> > rows;
		if constexpr (Matrix::Option_list::is_of_boundary_type){
			rows = build_longer_reduced_row_matrix<typename Matrix::Column_type>();
		} else {
			rows = build_longer_chain_row_matrix<typename Matrix::Column_type>();
		}
		unsigned int i = 0;
		for (auto& r : rows){
			if constexpr (Matrix::Option_list::has_removable_rows)
				if (i == 6) continue;
			test_column_equality<typename Matrix::Column_type>(r, get_ordered_row(m, i++));
		}
	}

	std::set<std::tuple<int,int,int>,BarComp> bars;
	//bars are not ordered the same for all matrices
	for (auto it = barcode.begin(); it != barcode.end(); ++it){
		bars.emplace(it->dim, it->birth, it->death);
	}
	auto it = bars.begin();
	BOOST_CHECK_EQUAL(std::get<0>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<2>(*it), -1);	//TODO: verify why this works...: it->death should be unsigned int, so double conversion
	++it;
	BOOST_CHECK_EQUAL(std::get<0>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 1);
	BOOST_CHECK_EQUAL(std::get<2>(*it), 3);
	++it;
	BOOST_CHECK_EQUAL(std::get<0>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 2);
	BOOST_CHECK_EQUAL(std::get<2>(*it), 4);
	++it;
	BOOST_CHECK_EQUAL(std::get<0>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 7);
	BOOST_CHECK_EQUAL(std::get<2>(*it), 8);
	++it;
	BOOST_CHECK_EQUAL(std::get<0>(*it), 1);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 5);
	BOOST_CHECK_EQUAL(std::get<2>(*it), 6);
	++it;
	BOOST_CHECK(it == bars.end());
}

// //*** swap/vine
// 	//base
// 	//boundary: does not update barcode
// 	void swap_columns(index columnIndex1, index columnIndex2);
// 	//base
// 	//boundary: does not update barcode
// 	void swap_rows(index rowIndex1, index rowIndex2);
// 	//base: no row access necessary (but then, has only effect with swaps), assumes the row is empty, just thought as an index cleanup
// 	void erase_row(index rowIndex);
// 	//*******************

template<class Matrix>
void test_base_swaps(){
	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
	Matrix m(columns, 5);

	test_content_equality(columns, m);

	m.swap_columns(3,5);
	if constexpr (Matrix::Option_list::column_indexation_type == Column_indexation_types::IDENTIFIER){
		m.swap_columns(3,1);
	} else {
		m.swap_columns(5,1);
		columns[3].swap(columns[5]);
		columns[5].swap(columns[1]);
	}

	test_content_equality(columns, m);

	m.swap_rows(3,5);
	m.swap_rows(5,1);

	if constexpr (Matrix::Option_list::is_z2){
		columns[4] = {2,5};
		columns[6] = {1,3,4};
		if constexpr (Matrix::Option_list::column_indexation_type == Column_indexation_types::IDENTIFIER){
			columns[3] = {0,5};
			columns[5] = {0,2};
		} else {
			columns[1] = {0,5};
		}
	} else {
		columns[4] = {{2,4},{5,1}};
		columns[6] = {{1,1},{3,4},{4,1}};
		if constexpr (Matrix::Option_list::column_indexation_type == Column_indexation_types::IDENTIFIER){
			columns[3] = {{0,1},{5,4}};
			columns[5] = {{0,1},{2,4}};
		} else {
			columns[1] = {{0,1},{5,4}};
		}
	}

	test_content_equality(columns, m);
}

// 	// //base
// 	// //boundary: does not update barcode
// 	// //id to pos
// 	// void swap_at_indices(index index1, index index2);
// 	// //*******************

// template<class Matrix>
// void test_base_index_swaps(){
// 	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
// 	Matrix m(columns, 5);

// 	test_content_equality(columns, m);

// 	m.swap_at_indices(3,5);
// 	m.swap_at_indices(5,1);

// 	columns[3].swap(columns[5]);
// 	columns[5].swap(columns[1]);
// 	if constexpr (Matrix::Option_list::is_z2){
// 		columns[1] = {0,5};
// 		columns[4] = {2,5};
// 		columns[6] = {1,3,4};
// 	} else {
// 		columns[1] = {{0,1},{5,4}};
// 		columns[4] = {{2,4},{5,1}};
// 		columns[6] = {{1,1},{3,4},{4,1}};
// 	}

// 	test_content_equality(columns, m);
// }

// template<class Matrix>
// void test_base_indexed_by_id_index_swaps(){
// 	auto columns = build_simple_boundary_matrix<typename Matrix::Column_type>();
// 	Matrix m(columns, 5);

// 	test_content_equality(columns, m);

// 	m.swap_at_indices(3,5);
// 	m.swap_at_indices(3,1);

// 	if constexpr (Matrix::Option_list::is_z2){
// 		columns[3] = {0,5};
// 		columns[4] = {2,5};
// 		columns[5] = {0,2};
// 		columns[6] = {1,3,4};
// 	} else {
// 		columns[3] = {{0,1},{5,4}};
// 		columns[4] = {{2,4},{5,1}};
// 		columns[5] = {{0,1},{2,4}};
// 		columns[6] = {{1,1},{3,4},{4,1}};
// 	}

// 	test_content_equality(columns, m);
// }

	// //ru: returns true if barcode was changed
	// //pos to id
	// bool vine_swap_with_z_eq_1_case(index index);								//by column position with ordered column container
	// //ru: returns true if barcode was changed
	// //pos to id
	// bool vine_swap(index index);												//by column position with ordered column container
	// //*******************

//assumes matrix was build with `build_longer_boundary_matrix` and was given the right comparision methods for non-barcode
template<class Matrix>
void test_vine_swap_with_position_index(Matrix& m){
	std::vector<witness_content<typename Matrix::Column_type> > columns;
	bool change;

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 8);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(5));
		BOOST_CHECK(m.is_zero_column(7));
		columns = build_longer_boundary_matrix<typename Matrix::Column_type>();
		columns[5].clear();
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 8);
		columns = build_longer_chain_matrix<typename Matrix::Column_type>();
	}
	test_content_equality(columns, m);

	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 3);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 5);
		BOOST_CHECK_EQUAL(it->death, 6);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 7);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK(it == barcode.end());
	}

	change = m.vine_swap(6);
	BOOST_CHECK(change);
	change = m.vine_swap(5);
	BOOST_CHECK(change);
	change = m.vine_swap(4);
	BOOST_CHECK(change);
	change = m.vine_swap(3);
	BOOST_CHECK(change);
	change = m.vine_swap(7);
	BOOST_CHECK(change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 8);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(3));
		BOOST_CHECK(m.is_zero_column(6));
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 7);
	}

	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 6);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 7);
		++it;
		BOOST_CHECK(it == barcode.end());
	}
	
	change = m.vine_swap(0);
	BOOST_CHECK(!change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 8);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(3));
		BOOST_CHECK(m.is_zero_column(6));
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 7);
	}
	
	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 6);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 7);
		++it;
		BOOST_CHECK(it == barcode.end());
	}

	change = m.vine_swap(4);
	BOOST_CHECK(change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 8);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(3));
		BOOST_CHECK(m.is_zero_column(6));
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 7);
	}
	
	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 6);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 7);
		++it;
		BOOST_CHECK(it == barcode.end());
	}

	change = m.vine_swap(5);
	BOOST_CHECK(!change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 8);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(3));
		BOOST_CHECK(m.is_zero_column(6));
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 7);
	}
	
	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 6);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 7);
		++it;
		BOOST_CHECK(it == barcode.end());
	}

	change = m.vine_swap(6);
	BOOST_CHECK(change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 8);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(3));
		BOOST_CHECK(m.is_zero_column(7));
		columns[3].clear();
		columns[4] = {0,2};
		columns[5] = {0,1};
		columns[6] = {0,3};
		columns[7].clear();
		columns[8] = {4,5,7};
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 6);
		columns[0] = {1};
		columns[1] = {0,1};
		columns[2] = {1,2};
		columns[3] = {0,7};
		columns[4] = {4};
		columns[5] = {4,5};
		columns[6] = {4,5,8};
		columns[7] = {3,4,5};
		columns[8] = {6};
	}
	test_content_equality(columns, m);
	
	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 7);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 6);
		++it;
		BOOST_CHECK(it == barcode.end());
	}
}

	// //chain: returns index which was not modified, ie new i+1
	// //id to pos
	// index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);	//by column id with potentielly unordered column container
	// //chain: returns index which was not modified, ie new i+1
	// //id to pos
	// index vine_swap(index columnIndex1, index columnIndex2);					//by column id with potentielly unordered column container
	// //*******************

//assumes matrix was build with `build_longer_boundary_matrix` and was given the right comparision methods for non-barcode
template<class Matrix>
void test_vine_swap_with_id_index(Matrix& m){
	std::vector<witness_content<typename Matrix::Column_type> > columns;
	unsigned int next;

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 8);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(5));
		BOOST_CHECK(m.is_zero_column(7));
		columns = build_longer_boundary_matrix<typename Matrix::Column_type>();
		columns[5].clear();
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 8);
		columns = build_longer_chain_matrix<typename Matrix::Column_type>();
	}
	test_content_equality(columns, m);
	
	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 3);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 5);
		BOOST_CHECK_EQUAL(it->death, 6);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 7);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK(it == barcode.end());
	}

	next = m.vine_swap(6, 7);
	BOOST_CHECK_EQUAL(next, 6);
	next = m.vine_swap(5, 7);
	BOOST_CHECK_EQUAL(next, 5);
	next = m.vine_swap(4, 7);
	BOOST_CHECK_EQUAL(next, 4);
	next = m.vine_swap(3, 7);
	BOOST_CHECK_EQUAL(next, 3);
	next = m.vine_swap(6, 8);
	BOOST_CHECK_EQUAL(next, 6);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(7));
		BOOST_CHECK(m.is_zero_column(5));
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 8);
	}
	
	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 6);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 7);
		++it;
		BOOST_CHECK(it == barcode.end());
	}

	next = m.vine_swap(0, 1);
	BOOST_CHECK_EQUAL(next, 1);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(7));
		BOOST_CHECK(m.is_zero_column(5));
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 8);
	}
	
	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 6);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 7);
		++it;
		BOOST_CHECK(it == barcode.end());
	}

	next = m.vine_swap(3, 4);
	BOOST_CHECK_EQUAL(next, 3);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(7));
		BOOST_CHECK(m.is_zero_column(5));
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 8);
	}
	
	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 6);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 7);
		++it;
		BOOST_CHECK(it == barcode.end());
	}

	next = m.vine_swap(3, 5);
	BOOST_CHECK_EQUAL(next, 5);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(7));
		BOOST_CHECK(m.is_zero_column(3));
	} else {
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 8);
	}

	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 6);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 7);
		++it;
		BOOST_CHECK(it == barcode.end());
	}

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		next = m.vine_swap(3, 8);		//use of simplex id
		BOOST_CHECK_EQUAL(next, 3);

		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 6);
		BOOST_CHECK(m.is_zero_column(0));
		BOOST_CHECK(m.is_zero_column(1));
		BOOST_CHECK(m.is_zero_column(2));
		BOOST_CHECK(m.is_zero_column(3));
		BOOST_CHECK(m.is_zero_column(7));

		columns[3].clear();
		columns[4] = {0,2};
		columns[5] = {0,1};
		columns[6] = {4,5,7};	//by id and not position
		columns[7].clear();
		columns[8] = {0,3};		//by id and not position
	} else {
		next = m.vine_swap(5, 8);		//use of internal chain id =/= initial simplex id
		BOOST_CHECK_EQUAL(next, 5);

		BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 5);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 3);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(m.get_column_with_pivot(8), 8);

		columns[0][0] = 1;
		columns[2][0] = 1;
		columns[3][0] = 4;
		columns[3].push_back(5);
		columns[4][0] = 4;
		columns[4].pop_back();
		columns[8][0] = 4;
		columns[8][1] = 5;
		columns[8].push_back(8);
	}
	test_content_equality(columns, m);
	
	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = m.get_current_barcode();
		auto it = barcode.begin();
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 0);
		BOOST_CHECK_EQUAL(it->death, static_cast<typename Matrix::pos_index>(-1));
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 1);
		BOOST_CHECK_EQUAL(it->death, 5);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 2);
		BOOST_CHECK_EQUAL(it->death, 4);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 1);
		BOOST_CHECK_EQUAL(it->birth, 7);
		BOOST_CHECK_EQUAL(it->death, 8);
		++it;
		BOOST_CHECK_EQUAL(it->dim, 0);
		BOOST_CHECK_EQUAL(it->birth, 3);
		BOOST_CHECK_EQUAL(it->death, 6);
		++it;
		BOOST_CHECK(it == barcode.end());
	}
}

// //*** Representative cycles
// 	//chain
// 	//ru
// 	//id to pos
// 	//pos to id
// 	void update_representative_cycles();
// 	//chain
// 	//ru
// 	//id to pos
// 	//pos to id
// 	const std::vector<cycle_type>& get_representative_cycles();
// 	//chain
// 	//ru
// 	//id to pos
// 	//pos to id
// 	const cycle_type& get_representative_cycle(const Bar& bar);

template<class Matrix>
void test_representative_cycles(Matrix& mb){
	mb.update_representative_cycles();

	const auto& cycles = mb.get_representative_cycles();
	BOOST_CHECK_EQUAL(cycles.size(), 5);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		std::vector<unsigned int> tmp;
		tmp.push_back(0);
		BOOST_CHECK(cycles[0] == tmp);
		tmp[0] = 1;
		BOOST_CHECK(cycles[1] == tmp);
		tmp[0] = 2;
		BOOST_CHECK(cycles[2] == tmp);
		tmp[0] = 3;
		tmp.push_back(4);
		tmp.push_back(5);
		BOOST_CHECK(cycles[3] == tmp);
		tmp.clear();
		tmp.push_back(7);
		BOOST_CHECK(cycles[4] == tmp);
	} else {
		std::vector<unsigned int> tmp;
		tmp.push_back(0);
		BOOST_CHECK(cycles[0] == tmp);
		tmp.push_back(1);
		BOOST_CHECK(cycles[1] == tmp);
		tmp[1] = 2;
		BOOST_CHECK(cycles[2] == tmp);
		tmp[0] = 3;
		tmp[1] = 4;
		tmp.push_back(5);
		BOOST_CHECK(cycles[3] == tmp);
		tmp[0] = 0;
		tmp[1] = 7;
		tmp.pop_back();
		BOOST_CHECK(cycles[4] == tmp);
	}

	if constexpr (Matrix::Option_list::has_column_pairings){
		const auto& barcode = mb.get_current_barcode();
		auto it = barcode.begin();
		for (auto& cycle : cycles){
			BOOST_CHECK(cycle == mb.get_representative_cycle(*it));
			++it;
		}
	}
}

#endif // PM_MATRIX_TESTS_H
