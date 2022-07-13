/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_CELL_LIST_COLUMN_H
#define Z2_CELL_LIST_COLUMN_H

#include "z2_reduced_cell_column_with_row.h"

namespace Gudhi {
namespace persistence_matrix {

struct z2_matrix_list_row_tag;
struct z2_matrix_list_column_tag;

using z2_base_hook_matrix_list_row = boost::intrusive::list_base_hook<
			boost::intrusive::tag < z2_matrix_list_row_tag >
		  , boost::intrusive::link_mode < boost::intrusive::auto_unlink >
	>;
using z2_base_hook_matrix_list_column = boost::intrusive::list_base_hook <
			boost::intrusive::tag < z2_matrix_list_column_tag >
		  , boost::intrusive::link_mode < boost::intrusive::safe_link >
		>;

struct Z2_list_cell : public Z2_row_cell, public z2_base_hook_matrix_list_row, public z2_base_hook_matrix_list_column
{
	Z2_list_cell(index columnIndex, index rowIndex)
		: Z2_row_cell(columnIndex, rowIndex){};
};

using Z2_list_column_type = boost::intrusive::list <
					Z2_list_cell
				  , boost::intrusive::constant_time_size<false>
				  , boost::intrusive::base_hook< z2_base_hook_matrix_list_column >  >;

using Z2_list_row_type = boost::intrusive::list <
					Z2_list_cell
				  , boost::intrusive::constant_time_size<false>
				  , boost::intrusive::base_hook< z2_base_hook_matrix_list_row >
				>;

template<class Master_matrix>
class Z2_reduced_cell_list_column_with_row
		: public Z2_reduced_cell_column_with_row<Z2_list_cell,
												 Z2_list_column_type,
												 Z2_list_row_type,
												 z2_base_hook_matrix_list_row,
												 typename Master_matrix::Column_pairing_option>
{
public:
	using Cell = Z2_list_cell;
	using Column_type = Z2_list_column_type;
	using Row_type = Z2_list_row_type;

	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	Z2_reduced_cell_list_column_with_row(matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
	template<class Chain_type>
	Z2_reduced_cell_list_column_with_row(index chainIndex, Chain_type& chain, dimension_type dimension, matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
	Z2_reduced_cell_list_column_with_row(const Z2_reduced_cell_list_column_with_row& other);

	void swap_independent_rows(index rowIndex);
	bool is_non_zero(index rowIndex);

	Z2_reduced_cell_list_column_with_row& operator+=(Z2_reduced_cell_list_column_with_row &column);
	template<class Friend_master_matrix>
	friend Z2_reduced_cell_list_column_with_row<Friend_master_matrix> operator+(
			Z2_reduced_cell_list_column_with_row<Friend_master_matrix> column1,
			Z2_reduced_cell_list_column_with_row<Friend_master_matrix> const& column2);

private:
	using RCC = Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,z2_base_hook_matrix_list_row,typename Master_matrix::Column_pairing_option>;

	matrix_type& matrix_;
	dictionnary_type& pivotToColumnIndex_;
};

template<class Master_matrix>
inline Z2_reduced_cell_list_column_with_row<Master_matrix>::Z2_reduced_cell_list_column_with_row(matrix_type& matrix, dictionnary_type& pivotToColumnIndex)
	: RCC(), matrix_(matrix), pivotToColumnIndex_(pivotToColumnIndex)
{}

template<class Master_matrix>
template<class Chain_type>
inline Z2_reduced_cell_list_column_with_row<Master_matrix>::Z2_reduced_cell_list_column_with_row(
		index chainIndex, Chain_type& chain, dimension_type dimension, matrix_type& matrix, dictionnary_type& pivotToColumnIndex)
	: RCC(chainIndex, chain, dimension), matrix_(matrix), pivotToColumnIndex_(pivotToColumnIndex)
{}

template<class Master_matrix>
inline Z2_reduced_cell_list_column_with_row<Master_matrix>::Z2_reduced_cell_list_column_with_row(
		const Z2_reduced_cell_list_column_with_row& other)
	: RCC(other), matrix_(other.matrix_), pivotToColumnIndex_(other.pivotToColumnIndex_)
{}

template<class Master_matrix>
inline void Z2_reduced_cell_list_column_with_row<Master_matrix>::swap_independent_rows(index rowIndex)
{
	std::swap(pivotToColumnIndex_.at(RCC::get_pivot()),
			  pivotToColumnIndex_.at(rowIndex));
	matrix_.at(pivotToColumnIndex_.at(RCC::get_pivot())).swap_rows(matrix_.at(pivotToColumnIndex_.at(rowIndex)));
}

template<class Master_matrix>
inline bool Z2_reduced_cell_list_column_with_row<Master_matrix>::is_non_zero(index rowIndex)
{
	for (Cell& cell : RCC::get_column())
		if (cell.get_row_index() == rowIndex) return true;

	return false;
}

template<class Master_matrix>
inline Z2_reduced_cell_list_column_with_row<Master_matrix> &Z2_reduced_cell_list_column_with_row<Master_matrix>::operator+=(Z2_reduced_cell_list_column_with_row &column)
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
		} else {
			if (it1->get_row_index() > it2->get_row_index()) {
				Cell * new_cell = new Cell(pos, it2->get_row_index());
				tc.insert(it1, *new_cell); //col link, in order
				matrix_.at(pivotToColumnIndex_.at(it2->get_row_index())).get_row().push_back(*new_cell);//row link,no order
				++it2;
			} else { //it1->key() == it2->key()
				typename Column_type::iterator tmp_it = it1;
				++it1;
				++it2;
				Cell* tmp_ptr = &(*tmp_it);
				tmp_it->z2_base_hook_matrix_list_row::unlink(); //unlink from row
				tc.erase(tmp_it); //remove from col
				delete tmp_ptr;
			}
		}
	}

	while (it2 != sc.end()) {//if it1 reached the end of its column, but not it2
		Cell * new_cell = new Cell(pos, it2->get_row_index());
		matrix_.at(pivotToColumnIndex_.at(it2->get_row_index())).get_row().push_back(*new_cell); //row links
		tc.insert(tc.end(), *new_cell);
		++it2;
	}

	if (!is_non_zero(RCC::get_lowest_simplex_index())){
		swap_independent_rows(column.get_pivot());
		RCC::swap_lowest_simplex_index(column);
	}

	return *this;
}

template<class Friend_master_matrix>
Z2_reduced_cell_list_column_with_row<Friend_master_matrix> operator+(
		Z2_reduced_cell_list_column_with_row<Friend_master_matrix> column1,
		Z2_reduced_cell_list_column_with_row<Friend_master_matrix> const& column2)
{
	column1 += column2;
	return column1;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_CELL_LIST_COLUMN_H
