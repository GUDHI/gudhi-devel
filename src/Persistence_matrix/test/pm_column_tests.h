/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_COLUMN_TESTS_H
#define PM_COLUMN_TESTS_H

#include <vector>
#include <set>

#include <boost/test/test_tools.hpp>

#include "pm_test_utilities.h"
#include "pm_column_tests_mastermatrix.h"

template<class Column>
using pool_type = Gudhi::persistence_matrix::New_cell_constructor<typename Column::Cell>;

//assumes column was not modified since construction, ie no duplicated / erased values in heap or vector column
template<class Column>
std::vector<column_content<Column> > get_ordered_column_contents(std::vector<Column>& matrix)
{
	std::vector<std::set<typename std::conditional<
					is_z2<Column>(),
					unsigned int,
					std::pair<unsigned int, typename Column::Field_element_type>
				>::type> > ordCol(matrix.size());
	for (unsigned int i = 0; i < matrix.size(); ++i){
		Column& col = matrix[i];
		for (auto& cell : col){
			if constexpr (is_z2<Column>()){
				ordCol[i].insert(cell.get_row_index());
			} else {
				ordCol[i].insert({cell.get_row_index(), cell.get_element()});
			}
		}
	}
	return ordCol;
}

//assumes column was not modified since construction, ie no duplicated / erased values in heap or vector column
template<class Column>
std::vector<column_content<Column> > get_ordered_rows(std::vector<Column>& matrix)
{
	std::vector<std::set<typename std::conditional<
					is_z2<Column>(),
					unsigned int,
					std::pair<unsigned int, typename Column::Field_element_type>
				>::type> > rows;
	for (Column& col : matrix){
		for (auto& cell : col){
			if (cell.get_row_index() >= rows.size()) rows.resize(cell.get_row_index() + 1);
			if constexpr (is_z2<Column>()){
				rows[cell.get_row_index()].insert(cell.get_column_index());
			} else {
				rows[cell.get_row_index()].insert({cell.get_column_index(), cell.get_element()});
			}
		}
	}
	return rows;
}

template<class Column>
std::vector<Column> build_column_matrix(pool_type<Column>* pool){
	std::vector<Column> matrix;

	if constexpr (is_z2<Column>()){
		using cont = std::vector<unsigned int>;
		matrix.emplace_back(cont{0,1,3,5}, 4, nullptr, pool);
		matrix.emplace_back(cont{0,1,2,5,6}, 4, nullptr, pool);
		matrix.emplace_back(cont{0,1,2,5,6}, 4, nullptr, pool);
		matrix.emplace_back(cont{}, 4, nullptr, pool);
		matrix.emplace_back(cont{0,1,3,5}, 4, nullptr, pool);
		matrix.emplace_back(matrix[1]);
	} else {
		using cont = std::vector<std::pair<unsigned int,typename Column::Field_element_type> >;
		matrix.emplace_back(cont{{0,1},{1,2},{3,3},{5,4}}, 4, &_g_operators, pool);
		matrix.emplace_back(cont{{0,4},{1,2},{2,1},{5,1},{6,1}}, 4, &_g_operators, pool);
		matrix.emplace_back(cont{{0,1},{1,3},{2,4},{5,4},{6,4}}, 4, &_g_operators, pool);
		matrix.emplace_back(cont{}, 4, &_g_operators, pool);
		matrix.emplace_back(cont{{0,1},{1,2},{3,3},{5,4}}, 4, &_g_operators, pool);
		matrix.emplace_back(matrix[1]);
	}

	return matrix;
}

template<class Column, class Rows>
std::vector<Column> build_column_matrix(Rows &rows, pool_type<Column>* pool){
	std::vector<Column> matrix;

	if constexpr (is_z2<Column>()){
		using cont = std::vector<unsigned int>;
		matrix.emplace_back(0, cont{0,1,3,5}, 4, &rows, nullptr, pool);
		matrix.emplace_back(1, cont{0,1,2,5,6}, 4, &rows, nullptr, pool);
		matrix.emplace_back(2, cont{0,1,2,5,6}, 4, &rows, nullptr, pool);
		matrix.emplace_back(3, cont{}, 4, &rows, nullptr, pool);
		matrix.emplace_back(4, cont{0,1,3,5}, 4, &rows, nullptr, pool);
		matrix.emplace_back(matrix[1], 5, &rows);
	} else {
		using cont = std::vector<std::pair<unsigned int,typename Column::Field_element_type> >;
		matrix.emplace_back(0, cont{{0,1},{1,2},{3,3},{5,4}}, 4, &rows, &_g_operators, pool);
		matrix.emplace_back(1, cont{{0,4},{1,2},{2,1},{5,1},{6,1}}, 4, &rows, &_g_operators, pool);
		matrix.emplace_back(2, cont{{0,1},{1,3},{2,4},{5,4},{6,4}}, 4, &rows, &_g_operators, pool);
		matrix.emplace_back(3, cont{}, 4, &rows, &_g_operators, pool);
		matrix.emplace_back(4, cont{{0,1},{1,2},{3,3},{5,4}}, 4, &rows, &_g_operators, pool);
		matrix.emplace_back(matrix[1], 5, &rows);
	}

	return matrix;
}

template<class Column, class Rows>
std::vector<Column> build_base_boundary_column_matrix(Rows &rows, pool_type<Column>* pool){
	std::vector<Column> matrix;

	if constexpr (is_z2<Column>()){
		using cont = std::vector<unsigned int>;
		matrix.emplace_back(0, cont{0,1,3,5}, &rows, nullptr, pool);
		matrix.emplace_back(1, cont{0,1,2,5,6}, &rows, nullptr, pool);
		matrix.emplace_back(2, cont{0,1,2,5,6}, &rows, nullptr, pool);
		matrix.emplace_back(3, cont{}, &rows, nullptr, pool);
		matrix.emplace_back(4, cont{0,1,3,5}, &rows, nullptr, pool);
		matrix.emplace_back(matrix[1], 5, &rows);
	} else {
		using cont = std::vector<std::pair<unsigned int,typename Column::Field_element_type> >;
		matrix.emplace_back(0, cont{{0,1},{1,2},{3,3},{5,4}}, &rows, &_g_operators, pool);
		matrix.emplace_back(1, cont{{0,4},{1,2},{2,1},{5,1},{6,1}}, &rows, &_g_operators, pool);
		matrix.emplace_back(2, cont{{0,1},{1,3},{2,4},{5,4},{6,4}}, &rows, &_g_operators, pool);
		matrix.emplace_back(3, cont{}, &rows, &_g_operators, pool);
		matrix.emplace_back(4, cont{{0,1},{1,2},{3,3},{5,4}}, &rows, &_g_operators, pool);
		matrix.emplace_back(matrix[1], 5, &rows);
	}

	return matrix;
}

template<class Column>
std::vector<Column> build_base_boundary_column_matrix(pool_type<Column>* pool){
	std::vector<Column> matrix;

	if constexpr (is_z2<Column>()){
		using cont = std::vector<unsigned int>;
		matrix.emplace_back(cont{0,1,3,5}, nullptr, pool);
		matrix.emplace_back(cont{0,1,2,5,6}, nullptr, pool);
		matrix.emplace_back(cont{0,1,2,5,6}, nullptr, pool);
		matrix.emplace_back(cont{}, nullptr, pool);
		matrix.emplace_back(cont{0,1,3,5}, nullptr, pool);
		matrix.emplace_back(matrix[1]);
	} else {
		using cont = std::vector<std::pair<unsigned int,typename Column::Field_element_type> >;
		matrix.emplace_back(cont{{0,1},{1,2},{3,3},{5,4}}, &_g_operators, pool);
		matrix.emplace_back(cont{{0,4},{1,2},{2,1},{5,1},{6,1}}, &_g_operators, pool);
		matrix.emplace_back(cont{{0,1},{1,3},{2,4},{5,4},{6,4}}, &_g_operators, pool);
		matrix.emplace_back(cont{}, &_g_operators, pool);
		matrix.emplace_back(cont{{0,1},{1,2},{3,3},{5,4}}, &_g_operators, pool);
		matrix.emplace_back(matrix[1]);
	}

	return matrix;
}

template<class Column>
std::vector<std::vector<typename std::conditional<
								  is_z2<Column>(),
								  unsigned int,
								  std::pair<unsigned int, typename Column::Field_element_type>
							   >::type> > build_rows()
{
	using cell_type = typename std::conditional<
								  is_z2<Column>(),
								  unsigned int,
								  std::pair<unsigned int, typename Column::Field_element_type>
							   >::type;
	using container_type = std::vector<cell_type>;
	std::vector<container_type> rows(7);

	if constexpr (is_z2<Column>()){
		rows[0] = {0,1,2,4,5};
		rows[1] = {0,1,2,4,5};
		rows[2] = {1,2,5};
		rows[3] = {0,4};
		rows[4] = {};
		rows[5] = {0,1,2,4,5};
		rows[6] = {1,2,5};
	} else {
		rows[0] = {{0,1},{1,4},{2,1},{4,1},{5,4}};
		rows[1] = {{0,2},{1,2},{2,3},{4,2},{5,2}};
		rows[2] = {{1,1},{2,4},{5,1}};
		rows[3] = {{0,3},{4,3}};
		rows[4] = {};
		rows[5] = {{0,4},{1,1},{2,4},{4,4},{5,1}};
		rows[6] = {{1,1},{2,4},{5,1}};
	}

	return rows;
}

template<class Column>
std::vector<std::vector<typename std::conditional<
								  is_z2<Column>(),
								  unsigned int,
								  std::pair<unsigned int, typename Column::Field_element_type>
							   >::type> > build_column_values()
{
	using cell_type = typename std::conditional<
								  is_z2<Column>(),
								  unsigned int,
								  std::pair<unsigned int, typename Column::Field_element_type>
							   >::type;
	using container_type = std::vector<cell_type>;
	std::vector<container_type> columns(6);

	if constexpr (is_z2<Column>()){
		columns[0] = {0,1,3,5};
		columns[1] = {0,1,2,5,6};
		columns[2] = {0,1,2,5,6};
		columns[3] = {};
		columns[4] = {0,1,3,5};
		columns[5] = {0,1,2,5,6};
	} else {
		columns[0] = {{0,1},{1,2},{3,3},{5,4}};
		columns[1] = {{0,4},{1,2},{2,1},{5,1},{6,1}};
		columns[2] = {{0,1},{1,3},{2,4},{5,4},{6,4}};
		columns[3] = {};
		columns[4] = {{0,1},{1,2},{3,3},{5,4}};
		columns[5] = {{0,4},{1,2},{2,1},{5,1},{6,1}};
	}

	return columns;
}

//**********constructors

// Intrusive_list_column();
// template<class Container_type>
// Intrusive_list_column(const Container_type& nonZeroChainRowIndices, dimension_type dimension);	//dimension gets ignored for base
// template<class Container_type, class Row_container_type>
// Intrusive_list_column(index columnIndex, const Container_type& nonZeroChainRowIndices, dimension_type dimension, Row_container_type &rowContainer);	//dimension gets ignored for base
// Intrusive_list_column(const Intrusive_list_column& column);
// template<class Row_container_type>
// Intrusive_list_column(const Intrusive_list_column& column, index columnIndex, Row_container_type &rowContainer);
// Intrusive_list_column(Intrusive_list_column&& column) noexcept;
// ~Intrusive_list_column();

// //Disabled with row access.
// Intrusive_list_column& operator=(const Intrusive_list_column& other);
// friend void swap(Intrusive_list_column& col1, Intrusive_list_column& col2){}

template<class Column>
void column_test_common_constructors(){
	using cell_type = typename std::conditional<
								  is_z2<Column>(),
								  unsigned int,
								  std::pair<unsigned int, typename Column::Field_element_type>
							   >::type;
	using container_type = std::vector<cell_type>;

	container_type cont1, cont2;
	typename Column::Field_operators *op = nullptr;
	pool_type<Column> pool;
	if constexpr (!is_z2<Column>()) op = &_g_operators;

	if constexpr (is_z2<Column>()){
		cont1 = {0,2,4};
		cont2 = {0,5,6};
	} else {
		cont1 = {{0, 1u},{2, 2u},{4, 3u}};
		cont2 = {{0, 1u},{5, 2u},{6, 3u}};
	}

	Column emptyCol(op, &pool);
	BOOST_CHECK_EQUAL(emptyCol.size(), 0);

	Column col(cont1, 2, op, &pool);
	BOOST_CHECK_EQUAL(col.size(), 3);

	std::vector<int> rows;	//type doesn't matter, as the row option is not enabled.
	Column rowCol(2, cont2, 2, &rows, op, &pool);
	BOOST_CHECK_EQUAL(rowCol.size(), 3);

	Column copyCol(col);
	BOOST_CHECK_EQUAL(copyCol.size(), 3);

	Column rowCopyCol(rowCol, 4, &rows);
	BOOST_CHECK_EQUAL(rowCopyCol.size(), 3);

	Column moveCol(std::move(col));
	BOOST_CHECK_EQUAL(moveCol.size(), 3);
	BOOST_CHECK_EQUAL(col.size(), 0);
	BOOST_CHECK(copyCol.get_content() == moveCol.get_content());

	Column assignCol = copyCol;
	BOOST_CHECK_EQUAL(assignCol.size(), 3);

	swap(copyCol, rowCol);
	BOOST_CHECK_EQUAL(copyCol.size(), 3);
	BOOST_CHECK_EQUAL(rowCol.size(), 3);
	BOOST_CHECK(copyCol.get_content() == rowCopyCol.get_content());
	BOOST_CHECK(rowCol.get_content() == assignCol.get_content());

	BOOST_CHECK_EQUAL(rows.size(), 0);
}

//**********methods

// std::vector<Field_element_type> get_content(int columnLength = -1) const;
// bool is_non_zero(index rowIndex) const;
// bool is_empty() const;
// std::size_t size() const;

// iterator begin() noexcept;
// const_iterator begin() const noexcept;
// iterator end() noexcept;
// const_iterator end() const noexcept;

template<class Column, typename cell_type>
void column_test_common_content_access(Column& col, 
									   const std::set<cell_type>& setcont, 
									   const std::vector<typename Column::Field_element_type>& veccont)
{
	BOOST_CHECK(get_column_content_via_iterators(col) == setcont);
	BOOST_CHECK(col.get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(col.size(), setcont.size());
	}

	if (veccont.empty()) BOOST_CHECK(col.is_empty());
	else BOOST_CHECK(!col.is_empty());

	for (unsigned int i = 0; i < veccont.size(); ++i){
		if (veccont[i] == 0u) BOOST_CHECK(!col.is_non_zero(i));
		else BOOST_CHECK(col.is_non_zero(i));
	}
}

//assumes that matrix was build with build_column_matrix and was not modified since.
template<class Column>
void column_test_common_z5_content_access(std::vector<Column> &matrix){
	std::set<std::pair<unsigned int,typename Column::Field_element_type> > setcont;
	std::vector<typename Column::Field_element_type> veccont;

	setcont = {{0,1},{1,2},{3,3},{5,4}};
	veccont = {1, 2, 0, 3, 0, 4};
	column_test_common_content_access(matrix[0], setcont, veccont);

	setcont = {{0,4},{1,2},{2,1},{5,1},{6,1}};
	veccont = {4, 2, 1, 0, 0, 1, 1};
	column_test_common_content_access(matrix[1], setcont, veccont);

	setcont = {{0,1},{1,3},{2,4},{5,4},{6,4}};
	veccont = {1, 3, 4, 0, 0, 4, 4};
	column_test_common_content_access(matrix[2], setcont, veccont);

	setcont = {};
	veccont = {};
	column_test_common_content_access(matrix[3], setcont, veccont);

	setcont = {{0,1},{1,2},{3,3},{5,4}};
	veccont = {1, 2, 0, 3, 0, 4};
	column_test_common_content_access(matrix[4], setcont, veccont);

	setcont = {{0,4},{1,2},{2,1},{5,1},{6,1}};
	veccont = {4, 2, 1, 0, 0, 1, 1};
	column_test_common_content_access(matrix[5], setcont, veccont);
}

//assumes that matrix was build with build_column_matrix and was not modified since.
template<class Column>
void column_test_common_z2_content_access(std::vector<Column> &matrix){
	std::set<unsigned int> setcont;
	std::vector<typename Column::Field_element_type> veccont;

	setcont = {0,1,3,5};
	veccont = {1, 1, 0, 1, 0, 1};
	column_test_common_content_access(matrix[0], setcont, veccont);

	setcont = {0,1,2,5,6};
	veccont = {1, 1, 1, 0, 0, 1, 1};
	column_test_common_content_access(matrix[1], setcont, veccont);

	setcont = {0,1,2,5,6};
	veccont = {1, 1, 1, 0, 0, 1, 1};
	column_test_common_content_access(matrix[2], setcont, veccont);

	setcont = {};
	veccont = {};
	column_test_common_content_access(matrix[3], setcont, veccont);

	setcont = {0,1,3,5};
	veccont = {1, 1, 0, 1, 0, 1};
	column_test_common_content_access(matrix[4], setcont, veccont);

	setcont = {0,1,2,5,6};
	veccont = {1, 1, 1, 0, 0, 1, 1};
	column_test_common_content_access(matrix[5], setcont, veccont);
}

// Intrusive_list_column& operator+=(Intrusive_list_column &column);	//for chain and vector
// friend Intrusive_list_column operator+(Intrusive_list_column column1, Intrusive_list_column& column2){}
// Intrusive_list_column& operator*=(unsigned int v);
// friend Intrusive_list_column operator*(Intrusive_list_column column, unsigned int const& v){}
// friend Intrusive_list_column operator*(unsigned int const& v, Intrusive_list_column column){}
// //this = v * this + column
// Intrusive_list_column& multiply_and_add(const Field_element_type& val, Intrusive_list_column& column);	//for chain and vector
// //this = this + column * v
// Intrusive_list_column& multiply_and_add(Intrusive_list_column& column, const Field_element_type& val);	//for chain and vector

// friend bool operator==(const Intrusive_list_column& c1, const Intrusive_list_column& c2){}
// friend bool operator<(const Intrusive_list_column& c1, const Intrusive_list_column& c2){}

//assumes that matrix was build with build_column_matrix and was not modified since.
template<class Column>
void column_test_common_z5_operators(std::vector<Column> &matrix){
	std::set<std::pair<unsigned int,typename Column::Field_element_type> > setcont;
	std::vector<typename Column::Field_element_type> veccont;

	matrix[0] += matrix[1];
	matrix[1] += matrix[2];

	setcont = {{1,4},{2,1},{3,3},{6,1}};
	veccont = {0, 4, 1, 3, 0, 0, 1};
	column_test_common_content_access(matrix[0], setcont, veccont);
	BOOST_CHECK(matrix[1] == matrix[1]);
	BOOST_CHECK(matrix[1] == matrix[3]);
	BOOST_CHECK(matrix[1] < matrix[0]);
	BOOST_CHECK(matrix[2] < matrix[0]);
	BOOST_CHECK(!(matrix[0] < matrix[0]));

	matrix[0] *= 4;
	matrix[1] *= 2;
	matrix[2] *= 1;

	setcont = {{1,1},{2,4},{3,2},{6,4}};
	veccont = {0, 1, 4, 2, 0, 0, 4};
	column_test_common_content_access(matrix[0], setcont, veccont);
	setcont = {{0,1},{1,3},{2,4},{5,4},{6,4}};
	veccont = {1, 3, 4, 0, 0, 4, 4};
	column_test_common_content_access(matrix[2], setcont, veccont);
	BOOST_CHECK(matrix[1] == matrix[1]);
	BOOST_CHECK(matrix[1] == matrix[3]);
	BOOST_CHECK(matrix[1] < matrix[0]);
	BOOST_CHECK(matrix[2] < matrix[0]);
	BOOST_CHECK(!(matrix[0] < matrix[0]));

	//this = v * this + column
	matrix[4].multiply_and_add(4, matrix[5]);
	setcont = {{0,3},{2,1},{3,2},{5,2},{6,1}};
	veccont = {3, 0, 1, 2, 0, 2, 1};
	column_test_common_content_access(matrix[4], setcont, veccont);
	//this = this + column * v
	matrix[5].multiply_and_add(matrix[3], 3);
	setcont = {{0,4},{1,2},{2,1},{5,1},{6,1}};
	veccont = {4, 2, 1, 0, 0, 1, 1};
	column_test_common_content_access(matrix[5], setcont, veccont);
	//this = this + column * v
	matrix[5].multiply_and_add(matrix[4], 3);
	setcont = {{0,3},{1,2},{2,4},{3,1},{5,2},{6,4}};
	veccont = {3, 2, 4, 1, 0, 2, 4};
	column_test_common_content_access(matrix[5], setcont, veccont);
	//this = v * this + column
	matrix[3].multiply_and_add(4, matrix[5]);
	setcont = {{0,3},{1,2},{2,4},{3,1},{5,2},{6,4}};
	veccont = {3, 2, 4, 1, 0, 2, 4};
	column_test_common_content_access(matrix[3], setcont, veccont);
	BOOST_CHECK(matrix[3] == matrix[5]);
	BOOST_CHECK(!(matrix[5] == matrix[4]));
	BOOST_CHECK(!(matrix[5] < matrix[3]));
	BOOST_CHECK(!(matrix[3] < matrix[5]));
	BOOST_CHECK(matrix[5] < matrix[4]);
	BOOST_CHECK(matrix[3] < matrix[4]);
}

//assumes that matrix was build with build_column_matrix and was not modified since.
template<class Column>
void column_test_common_z2_operators(std::vector<Column> &matrix){
	std::set<unsigned int> setcont;
	std::vector<typename Column::Field_element_type> veccont;

	matrix[0] += matrix[1];
	matrix[1] += matrix[2];

	setcont = {2,3,6};
	veccont = {0, 0, 1, 1, 0, 0, 1};
	column_test_common_content_access(matrix[0], setcont, veccont);
	BOOST_CHECK(matrix[1] == matrix[1]);
	BOOST_CHECK(matrix[1] == matrix[3]);
	BOOST_CHECK(matrix[1] < matrix[0]);
	BOOST_CHECK(matrix[2] < matrix[0]);
	BOOST_CHECK(!(matrix[0] < matrix[0]));

	matrix[0] *= 5;
	matrix[2] *= 1;

	setcont = {2,3,6};
	veccont = {0, 0, 1, 1, 0, 0, 1};
	column_test_common_content_access(matrix[0], setcont, veccont);
	setcont = {0,1,2,5,6};
	veccont = {1, 1, 1, 0, 0, 1, 1};
	column_test_common_content_access(matrix[2], setcont, veccont);
	BOOST_CHECK(!(matrix[0] == matrix[2]));
	BOOST_CHECK(matrix[2] == matrix[2]);
	BOOST_CHECK(matrix[1] < matrix[0]);
	BOOST_CHECK(matrix[2] < matrix[0]);
	BOOST_CHECK(!(matrix[0] < matrix[0]));

	//this = v * this + column
	matrix[4].multiply_and_add(3, matrix[5]);
	setcont = {2,3,6};
	veccont = {0, 0, 1, 1, 0, 0, 1};
	column_test_common_content_access(matrix[4], setcont, veccont);
	BOOST_CHECK(matrix[4] == matrix[0]);
	//this = this + column * v
	matrix[5].multiply_and_add(matrix[3], 3);
	setcont = {0,1,2,5,6};
	veccont = {1, 1, 1, 0, 0, 1, 1};
	column_test_common_content_access(matrix[5], setcont, veccont);
	BOOST_CHECK(matrix[5] == matrix[2]);
	//this = this + column * v
	matrix[5].multiply_and_add(matrix[4], 3);
	setcont = {0, 1,3,5};
	veccont = {1, 1, 0, 1, 0, 1, 0};
	column_test_common_content_access(matrix[5], setcont, veccont);
	//this = v * this + column
	matrix[3].multiply_and_add(3, matrix[5]);
	setcont = {0, 1,3,5};
	veccont = {1, 1, 0, 1, 0, 1, 0};
	column_test_common_content_access(matrix[3], setcont, veccont);
	BOOST_CHECK(matrix[3] == matrix[5]);
	BOOST_CHECK(!(matrix[5] < matrix[3]));
	BOOST_CHECK(!(matrix[3] < matrix[5]));
	BOOST_CHECK(matrix[5] < matrix[4]);
}

//**********row access

// template<class Container_type, class Row_container_type>
// Intrusive_list_column(index columnIndex, const Container_type& nonZeroChainRowIndices, dimension_type dimension, Row_container_type &rowContainer);	//dimension gets ignored for base
// template<class Row_container_type>
// Intrusive_list_column(const Intrusive_list_column& column, index columnIndex, Row_container_type &rowContainer);

//common: assumes that matrix was build with build_column_matrix(rows) and not modified since.
//base/boundary special: assumes that matrix was build with build_base_boundary_column_matrix(rows) and not modified since.
template<class Column, class row_container_type>
void column_test_row_access_constructors(std::vector<Column>& matrix, row_container_type& rows){
	test_matrix_equality<Column>(build_rows<Column>(), get_ordered_rows(matrix));
	test_matrix_equality<Column>(build_column_values<Column>(), get_ordered_column_contents(matrix));

	Column col(matrix[0], 6, &rows);
	for (auto& r : rows){
		if constexpr (Column::Master::Option_list::has_removable_rows){
			if (!r.second.empty()){
				auto& cell = *r.second.rbegin();
				if (cell.get_row_index() == 0 || cell.get_row_index() == 1 || cell.get_row_index() == 3 || cell.get_row_index() == 5){
					BOOST_CHECK_EQUAL(cell.get_column_index(), 6);
				}
			}
		} else {
			if (!r.empty()){
				auto& cell = *r.rbegin();
				if (cell.get_row_index() == 0 || cell.get_row_index() == 1 || cell.get_row_index() == 3 || cell.get_row_index() == 5){
					BOOST_CHECK_EQUAL(cell.get_column_index(), 6);
				}
			}
		}
	}

	Column moveCol(std::move(col));
	BOOST_CHECK_EQUAL(moveCol.size(), 4);
	BOOST_CHECK_EQUAL(col.size(), 0);
	BOOST_CHECK(matrix[0] == moveCol);
	for (auto& r : rows){
		if constexpr (Column::Master::Option_list::has_removable_rows){
			if (!r.second.empty()){
				auto& cell = *r.second.rbegin();
				if (cell.get_row_index() == 0 || cell.get_row_index() == 1 || cell.get_row_index() == 3 || cell.get_row_index() == 5){
					BOOST_CHECK_EQUAL(cell.get_column_index(), 6);
				}
			}
		} else {
			if (!r.empty()){
				auto& cell = *r.rbegin();
				if (cell.get_row_index() == 0 || cell.get_row_index() == 1 || cell.get_row_index() == 3 || cell.get_row_index() == 5){
					BOOST_CHECK_EQUAL(cell.get_column_index(), 6);
				}
			}
		}
	}

	swap(col, moveCol);
	BOOST_CHECK_EQUAL(moveCol.size(), 0);
	BOOST_CHECK_EQUAL(col.size(), 4);
	BOOST_CHECK(matrix[0] == col);
	for (auto& r : rows){
		if constexpr (Column::Master::Option_list::has_removable_rows){
			if (!r.second.empty()){
				auto& cell = *r.second.rbegin();
				if (cell.get_row_index() == 0 || cell.get_row_index() == 1 || cell.get_row_index() == 3 || cell.get_row_index() == 5){
					BOOST_CHECK_EQUAL(cell.get_column_index(), 6);
				}
			}
		} else {
			if (!r.empty()){
				auto& cell = *r.rbegin();
				if (cell.get_row_index() == 0 || cell.get_row_index() == 1 || cell.get_row_index() == 3 || cell.get_row_index() == 5){
					BOOST_CHECK_EQUAL(cell.get_column_index(), 6);
				}
			}
		}
	}
}

//**********constructors

// template<class Container_type>
// Intrusive_list_column(const Container_type& nonZeroRowIndices);	//has to be a boundary for boundary, has no sense for chain if dimension is needed
// template<class Container_type, class Row_container_type>
// Intrusive_list_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);	//has to be a boundary for boundary, has no sense for chain if dimension is needed

template<class Column>
void column_test_base_boundary_constructors(){
	using cell_type = typename std::conditional<
								  is_z2<Column>(),
								  unsigned int,
								  std::pair<unsigned int, typename Column::Field_element_type>
							   >::type;
	using container_type = std::vector<cell_type>;

	typename Column::Field_operators *op = nullptr;
	if constexpr (!is_z2<Column>()) op = &_g_operators;
	pool_type<Column> pool;

	container_type cont1, cont2;

	if constexpr (is_z2<Column>()){
		cont1 = {0,2,4};
		cont2 = {0,5,6};
	} else {
		cont1 = {{0, 1u},{2, 2u},{4, 3u}};
		cont2 = {{0, 1u},{5, 2u},{6, 3u}};
	}

	Column col(cont1, op, &pool);
	BOOST_CHECK_EQUAL(col.size(), 3);

	std::vector<int> rows;	//type doesn't matter, as the row option is not enabled.
	Column rowCol(2, cont2, &rows, op, &pool);
	BOOST_CHECK_EQUAL(rowCol.size(), 3);

	BOOST_CHECK(!(col == rowCol));
}

//**********methods

// //only for base and boundary
// template<class Map_type>
// void reorder(const Map_type& valueMap);	//used for lazy row swaps
// void clear();
// void clear(index rowIndex);

template<class Column>
void column_test_base_boundary_z5_methods(){
	std::vector<typename Column::Field_element_type> veccont;
	pool_type<Column> pool;

	Column col(std::vector<std::pair<unsigned int,typename Column::Field_element_type> >{{0,4},{1,2},{2,1},{5,1},{6,1}}, &_g_operators, &pool);
	veccont = {4, 2, 1, 0, 0, 1, 1};
	BOOST_CHECK(col.get_content(veccont.size()) == veccont);
	BOOST_CHECK_EQUAL(col.size(), 5);

	std::vector<unsigned int> permutation{0,5,2,3,1,4,6};
	col.reorder(permutation);
	veccont = {4, 0, 1, 0, 1, 2, 1};
	BOOST_CHECK(col.get_content(veccont.size()) == veccont);
	BOOST_CHECK_EQUAL(col.size(), 5);

	col.clear(2);
	col.clear(6);
	veccont = {4, 0, 0, 0, 1, 2, 0};
	BOOST_CHECK(col.get_content(veccont.size()) == veccont);
	BOOST_CHECK_EQUAL(col.size(), 3);

	col.clear();
	veccont = {};
	BOOST_CHECK(col.get_content(veccont.size()) == veccont);
	BOOST_CHECK_EQUAL(col.size(), 0);
}

template<class Column>
void column_test_base_boundary_z2_methods(){
	std::vector<typename Column::Field_element_type> veccont;
	pool_type<Column> pool;

	Column col(std::vector<unsigned int>{0,1,2,5,6}, nullptr, &pool);
	veccont = {1, 1, 1, 0, 0, 1, 1};
	BOOST_CHECK(col.get_content(veccont.size()) == veccont);
	BOOST_CHECK_EQUAL(col.size(), 5);

	std::vector<unsigned int> permutation{0,5,2,3,1,4,6};
	col.reorder(permutation);
	veccont = {1, 0, 1, 0, 1, 1, 1};
	BOOST_CHECK(col.get_content(veccont.size()) == veccont);
	BOOST_CHECK_EQUAL(col.size(), 5);

	col.clear(2);
	col.clear(6);
	veccont = {1, 0, 0, 0, 1, 1, 0};
	BOOST_CHECK(col.get_content(veccont.size()) == veccont);
	BOOST_CHECK_EQUAL(col.size(), 3);

	col.clear();
	veccont = {};
	BOOST_CHECK(col.get_content(veccont.size()) == veccont);
	BOOST_CHECK_EQUAL(col.size(), 0);
}

// template<class Cell_range>
// Intrusive_list_column& operator+=(const Cell_range& column);	//for base & boundary except vector

// //this = v * this + column
// template<class Cell_range>
// Intrusive_list_column& multiply_and_add(const Field_element_type& val, const Cell_range& column);	//for base & boundary except vector
// //this = this + column * v
// template<class Cell_range>
// Intrusive_list_column& multiply_and_add(const Cell_range& column, const Field_element_type& val);	//for base & boundary except vector

//assumes that matrix was build with build_column_matrix and was not modified since.
template<class Column>
void column_test_base_z5_operators(std::vector<Column> &matrix){
	using Cell = typename Column::Cell;
	std::set<Cell> setcont;
	std::vector<typename Column::Field_element_type> veccont;

	Cell cell(0);
	cell.set_element(4);
	setcont.insert(cell);
	cell = Cell(1);
	cell.set_element(2);
	setcont.insert(cell);
	cell = Cell(2);
	cell.set_element(1);
	setcont.insert(cell);
	cell = Cell(5);
	cell.set_element(1);
	setcont.insert(cell);
	cell = Cell(6);
	cell.set_element(1);
	setcont.insert(cell);
	matrix[0] += setcont;

	veccont = {0, 4, 1, 3, 0, 0, 1};
	BOOST_CHECK(matrix[0].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[0].size(), 4);
	}

	setcont.clear();
	cell = Cell(0);
	cell.set_element(1);
	setcont.insert(cell);
	cell = Cell(1);
	cell.set_element(3);
	setcont.insert(cell);
	cell = Cell(2);
	cell.set_element(4);
	setcont.insert(cell);
	cell = Cell(5);
	cell.set_element(4);
	setcont.insert(cell);
	cell = Cell(6);
	cell.set_element(4);
	setcont.insert(cell);
	matrix[1] += setcont;

	veccont = {};
	BOOST_CHECK(matrix[1].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[1].size(), 0);
	}

	//this = v * this + column
	matrix[4].multiply_and_add(4, setcont);
	veccont = {0, 1, 4, 2, 0, 0, 4};
	BOOST_CHECK(matrix[4].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[4].size(), 4);
	}
	//this = this + column * v
	setcont = {};
	matrix[5].multiply_and_add(setcont, 3);
	veccont = {4, 2, 1, 0, 0, 1, 1};
	BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[5].size(), 5);
	}
	//this = this + column * v
	setcont.clear();
	cell = Cell(0);
	cell.set_element(3);
	setcont.insert(cell);
	cell = Cell(2);
	cell.set_element(1);
	setcont.insert(cell);
	cell = Cell(3);
	cell.set_element(2);
	setcont.insert(cell);
	cell = Cell(5);
	cell.set_element(2);
	setcont.insert(cell);
	cell = Cell(6);
	cell.set_element(1);
	setcont.insert(cell);
	matrix[5].multiply_and_add(setcont, 3);
	veccont = {3, 2, 4, 1, 0, 2, 4};
	BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[5].size(), 6);
	}
	//this = v * this + column
	setcont.clear();
	cell = Cell(0);
	cell.set_element(3);
	setcont.insert(cell);
	cell = Cell(1);
	cell.set_element(2);
	setcont.insert(cell);
	cell = Cell(2);
	cell.set_element(4);
	setcont.insert(cell);
	cell = Cell(3);
	cell.set_element(1);
	setcont.insert(cell);
	cell = Cell(5);
	cell.set_element(2);
	setcont.insert(cell);
	cell = Cell(6);
	cell.set_element(4);
	setcont.insert(cell);
	matrix[3].multiply_and_add(4, setcont);
	veccont = {3, 2, 4, 1, 0, 2, 4};
	BOOST_CHECK(matrix[3].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[3].size(), 6);
	}
}

//assumes that matrix was build with build_column_matrix and was not modified since.
template<class Column>
void column_test_base_z2_operators(std::vector<Column> &matrix){
	std::set<typename Column::Cell> setcont;
	std::vector<bool> veccont;

	setcont = {0,1,2,5,6};
	matrix[0] += setcont;

	veccont = {0, 0, 1, 1, 0, 0, 1};
	BOOST_CHECK(matrix[0].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[0].size(), 3);
	}

	setcont = {0,1,2,5,6};
	matrix[1] += setcont;

	veccont = {};
	BOOST_CHECK(matrix[1].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[1].size(), 0);
	}

	//this = v * this + column
	matrix[4].multiply_and_add(1, setcont);
	veccont = {0, 0, 1, 1, 0, 0, 1};
	BOOST_CHECK(matrix[4].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[4].size(), 3);
	}
	//this = this + column * v
	setcont = {};
	matrix[5].multiply_and_add(setcont, 1);
	veccont = {1, 1, 1, 0, 0, 1, 1};
	BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[5].size(), 5);
	}
	//this = this + column * v
	setcont = {2,3,6};
	matrix[5].multiply_and_add(setcont, 1);
	veccont = {1, 1, 0, 1, 0, 1, 0};
	BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[5].size(), 4);
	}
	//this = v * this + column
	setcont = {0,1,3,5};
	matrix[3].multiply_and_add(0, setcont);
	veccont = {1, 1, 0, 1, 0, 1, 0};
	BOOST_CHECK(matrix[3].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[3].size(), 4);
	}
}

//**********row access

// template<class Container_type, class Row_container_type>
// Intrusive_list_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);	//has to be a boundary for boundary, has no sense for chain if dimension is needed

//**********constructors

// template<class Container_type>
// Intrusive_list_column(const Container_type& nonZeroRowIndices);	//has to be a boundary for boundary, has no sense for chain if dimension is needed
// template<class Container_type, class Row_container_type>
// Intrusive_list_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);	//has to be a boundary for boundary, has no sense for chain if dimension is needed

//**********methods

// //only for base and boundary
// template<class Map_type>
// void reorder(const Map_type& valueMap);	//used for lazy row swaps
// void clear();
// void clear(index rowIndex);

// //only for chain and boundary
// int get_pivot() const;
// Field_element_type get_pivot_value() const;
// dimension_type get_dimension() const { return dim_; }

template<class Column>
void column_test_boundary_chain_methods(std::vector<Column>& matrix){
	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 5);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 6);
	BOOST_CHECK_EQUAL(matrix[2].get_pivot(), 6);
	BOOST_CHECK_EQUAL(matrix[3].get_pivot(), -1);
	BOOST_CHECK_EQUAL(matrix[4].get_pivot(), 5);
	BOOST_CHECK_EQUAL(matrix[5].get_pivot(), 6);

	if constexpr (!is_z2<Column>()){
		BOOST_CHECK_EQUAL(matrix[0].get_pivot_value(), 4u);
		BOOST_CHECK_EQUAL(matrix[1].get_pivot_value(), 1u);
		BOOST_CHECK_EQUAL(matrix[2].get_pivot_value(), 4u);
		BOOST_CHECK_EQUAL(matrix[3].get_pivot_value(), 0u);
		BOOST_CHECK_EQUAL(matrix[4].get_pivot_value(), 4u);
		BOOST_CHECK_EQUAL(matrix[5].get_pivot_value(), 1u);
	}

	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[2].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[3].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[4].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[5].get_dimension(), 4);
}

template<class Column>
void column_test_boundary_methods(std::vector<Column>& matrix){
	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 3);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[2].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[3].get_dimension(), 0);
	BOOST_CHECK_EQUAL(matrix[4].get_dimension(), 3);
	BOOST_CHECK_EQUAL(matrix[5].get_dimension(), 4);
}

// Intrusive_list_column& operator+=(const Intrusive_list_column& column);	//for base & boundary except vector

// //this = v * this + column
// Intrusive_list_column& multiply_and_add(const Field_element_type& val, const Intrusive_list_column& column);	//for base & boundary except vector
// //this = this + column * v
// Intrusive_list_column& multiply_and_add(const Intrusive_list_column& column, const Field_element_type& val);	//for base & boundary except vector

//assumes that matrix was build with build_column_matrix and was not modified since.
template<class Column>
void column_test_boundary_z5_operators(std::vector<Column> &matrix){
	using cont = std::vector<std::pair<unsigned int,typename Column::Field_element_type> >;

	std::vector<typename Column::Field_element_type> veccont;
	pool_type<Column> pool;

	const Column col0(cont{{0,4},{1,2},{2,1},{5,1},{6,1}}, &_g_operators, &pool);
	matrix[0] += col0;

	veccont = {0, 4, 1, 3, 0, 0, 1};
	BOOST_CHECK(matrix[0].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[0].size(), 4);
	}

	const Column col1(cont{{0,1},{1,3},{2,4},{5,4},{6,4}}, &_g_operators, &pool);
	matrix[1] += col1;

	veccont = {};
	BOOST_CHECK(matrix[1].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[1].size(), 0);
	}

	//this = v * this + column
	matrix[4].multiply_and_add(4, col1);
	veccont = {0, 1, 4, 2, 0, 0, 4};
	BOOST_CHECK(matrix[4].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[4].size(), 4);
	}
	//this = this + column * v
	const Column col2(cont{}, &_g_operators, &pool);
	matrix[5].multiply_and_add(col2, 3);
	veccont = {4, 2, 1, 0, 0, 1, 1};
	BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
	BOOST_CHECK_EQUAL(matrix[5].size(), 5);
	//this = this + column * v
	const Column col3(cont{{0,3},{2,1},{3,2},{5,2},{6,1}}, &_g_operators, &pool);
	matrix[5].multiply_and_add(col3, 3);
	veccont = {3, 2, 4, 1, 0, 2, 4};
	BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[5].size(), 6);
	}
	//this = v * this + column
	const Column col4(cont{{0,3},{1,2},{2,4},{3,1},{5,2},{6,4}}, &_g_operators, &pool);
	matrix[3].multiply_and_add(4, col4);
	veccont = {3, 2, 4, 1, 0, 2, 4};
	BOOST_CHECK(matrix[3].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[3].size(), 6);
	}
}

//assumes that matrix was build with build_column_matrix and was not modified since.
template<class Column>
void column_test_boundary_z2_operators(std::vector<Column> &matrix){
	using cont = std::vector<unsigned int>;

	std::vector<bool> veccont;
	pool_type<Column> pool;

	const Column col0(cont{0,1,2,5,6}, nullptr, &pool);
	matrix[0] += col0;

	veccont = {0, 0, 1, 1, 0, 0, 1};
	BOOST_CHECK(matrix[0].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[0].size(), 3);
	}

	const Column col1(cont{0,1,2,5,6}, nullptr, &pool);
	matrix[1] += col1;

	veccont = {};
	BOOST_CHECK(matrix[1].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[1].size(), 0);
	}

	//this = v * this + column
	matrix[4].multiply_and_add(3, col1);
	veccont = {0, 0, 1, 1, 0, 0, 1};
	BOOST_CHECK(matrix[4].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[4].size(), 3);
	}
	//this = this + column * v
	const Column col2(cont{}, nullptr, &pool);
	matrix[5].multiply_and_add(col2, 3);
	veccont = {1, 1, 1, 0, 0, 1, 1};
	BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[5].size(), 5);
	}
	//this = this + column * v
	const Column col3(cont{2,3,6}, nullptr, &pool);
	matrix[5].multiply_and_add(col3, 3);
	veccont = {1, 1, 0, 1, 0, 1, 0};
	BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[5].size(), 4);
	}
	//this = v * this + column
	const Column col4(cont{0,1,3,5}, nullptr, &pool);
	matrix[3].multiply_and_add(4, col4);
	veccont = {1, 1, 0, 1, 0, 1, 0};
	BOOST_CHECK(matrix[3].get_content(veccont.size()) == veccont);
	if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP){
		BOOST_CHECK_EQUAL(matrix[3].size(), 4);
	}
}

//**********row access

// template<class Container_type, class Row_container_type>
// Intrusive_list_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);	//has to be a boundary for boundary, has no sense for chain if dimension is needed

//**********constructors

//**********methods

// //only for chain and boundary
// int get_pivot() const;
// Field_element_type get_pivot_value() const;
// dimension_type get_dimension() const { return dim_; }

// int get_paired_chain_index() const { return pairedColumn_; }
// bool is_paired() const { return pairedColumn_ != -1; }
// void assign_paired_chain(index other_col){ pairedColumn_ = other_col; }
// void unassign_paired_chain() { pairedColumn_ = -1; };

template<class Column>
void column_test_chain_methods(){
	typename Column::Field_operators *op = nullptr;
	if constexpr (!is_z2<Column>()) op = &_g_operators;
	pool_type<Column> pool;

	Column col(op, &pool);

	BOOST_CHECK(!col.is_paired());
	BOOST_CHECK(col.get_paired_chain_index() == static_cast<typename Column::index>(-1));

	col.unassign_paired_chain();
	BOOST_CHECK(!col.is_paired());
	BOOST_CHECK(col.get_paired_chain_index() == static_cast<typename Column::index>(-1));

	col.assign_paired_chain(2);
	BOOST_CHECK(col.is_paired());
	BOOST_CHECK(col.get_paired_chain_index() == 2);

	col.unassign_paired_chain();
	BOOST_CHECK(!col.is_paired());
	BOOST_CHECK(col.get_paired_chain_index() == static_cast<typename Column::index>(-1));
}

#endif // PM_COLUMN_TESTS_H
