/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CELL_LIST_COLUMN_H
#define CELL_LIST_COLUMN_H

#include "reduced_cell_column_with_row.h"

namespace Gudhi {
namespace persistence_matrix {

struct matrix_list_row_tag;
struct matrix_list_column_tag;

using base_hook_matrix_list_row = boost::intrusive::list_base_hook<
			boost::intrusive::tag < matrix_list_row_tag >
		  , boost::intrusive::link_mode < boost::intrusive::auto_unlink >
	>;
using base_hook_matrix_list_column = boost::intrusive::list_base_hook <
			boost::intrusive::tag < matrix_list_column_tag >
		  , boost::intrusive::link_mode < boost::intrusive::safe_link >
		>;

template<class Field_element_type>
struct List_cell : public Row_cell<Field_element_type>, public base_hook_matrix_list_row, public base_hook_matrix_list_column
{
	List_cell(unsigned int element, index columnIndex, index rowIndex)
		: Row_cell<Field_element_type>(element, columnIndex, rowIndex){};
};

template<class Field_element_type>
using List_column_type = boost::intrusive::list <
					List_cell<Field_element_type>
				  , boost::intrusive::constant_time_size<false>
				  , boost::intrusive::base_hook< base_hook_matrix_list_column >  >;

template<class Field_element_type>
using List_row_type = boost::intrusive::list <
					List_cell<Field_element_type>
				  , boost::intrusive::constant_time_size<false>
				  , boost::intrusive::base_hook< base_hook_matrix_list_row >
				>;

template<class Master_matrix>
class Reduced_cell_list_column_with_row
		: public Reduced_cell_column_with_row<List_cell<typename Master_matrix::Field_type>,
											  List_column_type<typename Master_matrix::Field_type>,
											  List_row_type<typename Master_matrix::Field_type>,
											  base_hook_matrix_list_row,
											  typename Master_matrix::Field_type,
											  typename Master_matrix::Column_pairing_option>
{
public:
	using Field_element_type = typename Master_matrix::Field_type;
	using Cell = List_cell<Field_element_type>;
	using Column_type = List_column_type<Field_element_type>;
	using Row_type = List_row_type<Field_element_type>;

	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	Reduced_cell_list_column_with_row(matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
	template<class Chain_type>
	Reduced_cell_list_column_with_row(index chainIndex, Chain_type& chain, dimension_type dimension, matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
	Reduced_cell_list_column_with_row(const Reduced_cell_list_column_with_row& other);

	void swap_independent_rows(index rowIndex);
	bool is_non_zero(index rowIndex);

	Reduced_cell_list_column_with_row& operator+=(Reduced_cell_list_column_with_row &column);
	template<class Friend_master_matrix>
	friend Reduced_cell_list_column_with_row<Friend_master_matrix> operator+(
			Reduced_cell_list_column_with_row<Friend_master_matrix> column1,
			Reduced_cell_list_column_with_row<Friend_master_matrix> const& column2);
	Reduced_cell_list_column_with_row& operator*=(unsigned int v);
	template<class Friend_master_matrix>
	friend Reduced_cell_list_column_with_row<Friend_master_matrix> operator*(
			Reduced_cell_list_column_with_row<Friend_master_matrix> column,
			unsigned int const& v);
	template<class Friend_master_matrix>
	friend Reduced_cell_list_column_with_row<Friend_master_matrix> operator*(
			unsigned int const& v,
			Reduced_cell_list_column_with_row<Friend_master_matrix> column);

private:
	using RCC = Reduced_cell_column_with_row<Cell,Column_type,Row_type,base_hook_matrix_list_row,Field_element_type,typename Master_matrix::Column_pairing_option>;

	matrix_type& matrix_;
	dictionnary_type& pivotToColumnIndex_;
};

template<class Master_matrix>
inline Reduced_cell_list_column_with_row<Master_matrix>::Reduced_cell_list_column_with_row(matrix_type& matrix, dictionnary_type& pivotToColumnIndex)
	: RCC(), matrix_(matrix), pivotToColumnIndex_(pivotToColumnIndex)
{}

template<class Master_matrix>
template<class Chain_type>
inline Reduced_cell_list_column_with_row<Master_matrix>::Reduced_cell_list_column_with_row(
		index chainIndex, Chain_type& chain, dimension_type dimension, matrix_type& matrix, dictionnary_type& pivotToColumnIndex)
	: RCC(chainIndex, chain, dimension), matrix_(matrix), pivotToColumnIndex_(pivotToColumnIndex)
{}

template<class Master_matrix>
inline Reduced_cell_list_column_with_row<Master_matrix>::Reduced_cell_list_column_with_row(
		const Reduced_cell_list_column_with_row& other)
	: RCC(other), matrix_(other.matrix_), pivotToColumnIndex_(other.pivotToColumnIndex_)
{}

template<class Master_matrix>
inline void Reduced_cell_list_column_with_row<Master_matrix>::swap_independent_rows(index rowIndex)
{
	std::swap(pivotToColumnIndex_.at(RCC::get_pivot()),
			  pivotToColumnIndex_.at(rowIndex));
	matrix_.at(pivotToColumnIndex_.at(RCC::get_pivot())).swap_rows(matrix_.at(pivotToColumnIndex_.at(rowIndex)));
}

template<class Master_matrix>
inline bool Reduced_cell_list_column_with_row<Master_matrix>::is_non_zero(index rowIndex)
{
	for (Cell& cell : RCC::get_column())
		if (cell.get_row_index() == rowIndex) return true;

	return false;
}

template<class Master_matrix>
inline Reduced_cell_list_column_with_row<Master_matrix> &Reduced_cell_list_column_with_row<Master_matrix>::operator+=(Reduced_cell_list_column_with_row &column)
{
	Column_type& tc = RCC::get_column();
	Column_type& sc = column.get_column();
	index pos = pivotToColumnIndex_.at(RCC::get_pivot());

	auto it1 = tc.begin();
	auto it2 = sc.begin();
	while (it1 != tc.end() && it2 != sc.end())
	{
		if (it1->get_row_index() < it2->get_row_index()) {
			++it1;
		} else if (it1->get_row_index() > it2->get_row_index()) {
			Cell* new_cell = new Cell(it2->get_element(), pos, it2->get_row_index());
			tc.insert(it1, *new_cell);
			matrix_.at(pivotToColumnIndex_.at(it2->get_row_index())).get_row().push_back(*new_cell);
			++it2;
		} else {
			it1->get_element() += it2->get_element();
			if (it1->get_element() == 0u){
				typename Column_type::iterator tmp_it = it1;
				++it1;
				++it2;
				Cell* tmp_ptr = &(*tmp_it);
				tmp_it->base_hook_matrix_list_row::unlink();
				tc.erase(tmp_it);
				delete tmp_ptr;
			}
		}
	}

	while (it2 != sc.end()) {
		Cell* new_cell = new Cell(it2->get_element(), pos, it2->get_row_index());
		matrix_.at(pivotToColumnIndex_.at(it2->get_row_index())).get_row().push_back(*new_cell);
		tc.insert(tc.end(), *new_cell);
		++it2;
	}

	if (!is_non_zero(RCC::get_lowest_simplex_index())){
		swap_independent_rows(column.get_pivot());
		RCC::swap_lowest_simplex_index(column);
	}

	return *this;
}

template<class Master_matrix>
inline Reduced_cell_list_column_with_row<Master_matrix> &Reduced_cell_list_column_with_row<Master_matrix>::operator*=(unsigned int v)
{
	v %= Field_element_type::get_characteristic();

	if (v == 0) {
		auto it = RCC::get_column().begin();
		while (it != RCC::get_column().end()){
			typename Column_type::iterator tmp_it = it;
			++it;
			Cell* tmp_ptr = &(*tmp_it);
			tmp_it->base_hook_matrix_list_row::unlink();
			RCC::get_column().erase(tmp_it);
			delete tmp_ptr;
		}
		return *this;
	}

	for (Cell& cell : RCC::get_column()){
		cell.get_element() *= v;
	}

	return *this;
}

template<class Friend_master_matrix>
Reduced_cell_list_column_with_row<Friend_master_matrix> operator+(
		Reduced_cell_list_column_with_row<Friend_master_matrix> column1,
		Reduced_cell_list_column_with_row<Friend_master_matrix> const& column2)
{
	column1 += column2;
	return column1;
}

template<class Friend_master_matrix>
Reduced_cell_list_column_with_row<Friend_master_matrix> operator*(
		Reduced_cell_list_column_with_row<Friend_master_matrix> column,
		unsigned int const& v)
{
	column *= v;
	return column;
}

template<class Friend_master_matrix>
Reduced_cell_list_column_with_row<Friend_master_matrix> operator*(
		unsigned int const& v,
		Reduced_cell_list_column_with_row<Friend_master_matrix> column)
{
	column *= v;
	return column;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CELL_LIST_COLUMN_H
