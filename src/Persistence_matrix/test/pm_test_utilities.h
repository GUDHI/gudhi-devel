/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_TEST_UTILITIES_H
#define PM_TEST_UTILITIES_H

#include <type_traits>	//std::is_same_v
#include <set>

#include <boost/test/test_tools.hpp>
#include <gudhi/persistence_matrix_options.h>
#include <gudhi/Persistence_matrix/columns/heap_column.h>
#include <gudhi/matrix.h>
#include <gudhi/Fields/Zp_field_operators.h>

using Gudhi::persistence_matrix::Heap_column;
using Gudhi::persistence_matrix::Matrix;
using Gudhi::persistence_matrix::Column_indexation_types;
using Zp = Gudhi::persistence_fields::Zp_field_operators<> ;

inline Zp _g_operators = Zp(5);

template<class Column>
constexpr bool is_z2(){
	return std::is_same_v<typename Column::Field_element_type, bool>;
}

template<class Column>
constexpr bool has_row_access(){
	return Column::Master::Option_list::has_row_access;
}

template<class Matrix>
constexpr bool is_Boundary(){
	return Matrix::isNonBasic && Matrix::Option_list::is_of_boundary_type && !Matrix::Option_list::has_vine_update & !Matrix::Option_list::can_retrieve_representative_cycles;
}

template<class Matrix>
constexpr bool is_RU(){
	return Matrix::isNonBasic && Matrix::Option_list::is_of_boundary_type && (Matrix::Option_list::has_vine_update || Matrix::Option_list::can_retrieve_representative_cycles);
}

template<class Matrix>
constexpr bool is_Chain(){
	return Matrix::isNonBasic && !Matrix::Option_list::is_of_boundary_type;
}

template<class Matrix>
constexpr bool is_indexed_by_position(){
	return Matrix::Option_list::column_indexation_type == Column_indexation_types::POSITION || 
			(!is_Chain<Matrix>() && Matrix::Option_list::column_indexation_type == Column_indexation_types::CONTAINER);
}

template<class Column>
using cell_rep_type = typename std::conditional<
							is_z2<Column>(),
							unsigned int,
							std::pair<unsigned int, typename Column::Field_element_type>
						>::type;

template<class Column>
using column_content = std::set<cell_rep_type<Column> >;
template<class Column>
using witness_content = std::vector<cell_rep_type<Column> >;

//for vector, assumes no cell was removed via clear(index)
template<class Column>
column_content<Column> get_column_content_via_iterators(const Column& col)
{
	column_content<Column> cells;

	for (const auto& c : col){
		if constexpr (is_z2<Column>()){
			auto p = cells.insert(c.get_row_index());	//for vector, assumes no clear(idx) was called
		} else {
			auto p = cells.insert({c.get_row_index(), c.get_element()});	//for vector, assumes no clear(idx) was called
		}
	}

	return cells;
}

//assumes no cell was removed via clear(index)
template<class Matrix>
column_content<Heap_column<Matrix> > get_column_content_via_iterators(const Heap_column<Matrix>& col)
{
	column_content<Heap_column<Matrix> > cells;
	std::vector<typename Heap_column<Matrix>::Field_element_type> cont;

	for (const auto& c : col){
		if constexpr (is_z2<Heap_column<Matrix> >()){
			auto p = cells.insert(c.get_row_index());
			if (!p.second){	//possible in heap
				cells.erase(p.first);
			}
		} else {
			if (cont.size() <= c.get_row_index()) cont.resize(c.get_row_index() + 1, 0u);
			cont[c.get_row_index()] = _g_operators.add(cont[c.get_row_index()], c.get_element());
		}
	}
	if constexpr (!is_z2<Heap_column<Matrix> >()){
		for (unsigned int i = 0; i < cont.size(); ++i){
			if (cont[i] != 0u) cells.insert({i, cont[i]});
		}
	}

	return cells;
}

//for vector, assumes no cell was removed via clear(index)
template<class Matrix>
column_content<typename Matrix::Column_type> get_ordered_row(Matrix& m, unsigned int rowIndex)	//base and base comp cannot call get_row as const...
{
	column_content<typename Matrix::Column_type> orderedRows;
	for (const auto& cell : m.get_row(rowIndex)){
		if constexpr (is_z2<typename Matrix::Column_type>()){
			orderedRows.insert(cell.get_column_index());
		} else {
			orderedRows.insert({cell.get_column_index(), cell.get_element()});
		}
	}
	return orderedRows;
}

template<class Column>
void test_column_equality(const witness_content<Column>& witness, const column_content<Column>& ord)
{
	BOOST_CHECK_EQUAL(ord.size(), witness.size());
	auto itW = witness.begin();
	auto itO = ord.begin();
	while (itW != witness.end()){
		if constexpr (is_z2<Column>()){
			BOOST_CHECK_EQUAL(*itW, *itO);
		} else {
			BOOST_CHECK_EQUAL(itW->first, itO->first);
			BOOST_CHECK_EQUAL(itW->second, itO->second);
		}
		++itW; ++itO;
	}
}

template<class Column>
void test_matrix_equality(const std::vector<witness_content<Column> >& witness, 
						  const std::vector<column_content<Column> >& ord)
{
	BOOST_CHECK_EQUAL(witness.size(), ord.size());
	for (unsigned int i = 0; i < witness.size(); ++i){
		test_column_equality<Column>(witness[i], ord[i]);
	}
}

template<class Matrix>
void test_content_equality(const std::vector<witness_content<typename Matrix::Column_type> >& witness, 
						   Matrix& m)	//base and base comp cannot call get_column as const...
{
	unsigned int i = 0;
	for (auto& b : witness){
		const auto& col = m.get_column(i++);	//to force the const version
		test_column_equality<typename Matrix::Column_type>(b, get_column_content_via_iterators(col));
	}
}



#endif // PM_TEST_UTILITIES_H
