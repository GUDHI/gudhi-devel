/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_CELL_SET_COLUMN_H
#define Z2_CELL_SET_COLUMN_H

#include "z2_reduced_cell_column_with_row.h"

namespace Gudhi {
namespace persistence_matrix {

struct z2_matrix_set_row_tag;
struct z2_matrix_set_column_tag;

using z2_base_hook_matrix_set_row = boost::intrusive::list_base_hook<
			boost::intrusive::tag < z2_matrix_set_row_tag >
		  , boost::intrusive::link_mode < boost::intrusive::auto_unlink >
	>;
using z2_base_hook_matrix_set_column = boost::intrusive::set_base_hook <
			boost::intrusive::tag < z2_matrix_set_column_tag >
		  , boost::intrusive::link_mode < boost::intrusive::safe_link >
		>;

struct Z2_set_cell : public Z2_row_cell, public z2_base_hook_matrix_set_row, public z2_base_hook_matrix_set_column
{
	Z2_set_cell(index columnIndex, index rowIndex)
		: Z2_row_cell(columnIndex, rowIndex){};
};

using Z2_set_column_type = boost::intrusive::set <
					Z2_set_cell
				  , boost::intrusive::constant_time_size<false>
				  , boost::intrusive::base_hook< z2_base_hook_matrix_set_column >  >;

using Z2_set_row_type = boost::intrusive::list <
					Z2_set_cell
				  , boost::intrusive::constant_time_size<false>
				  , boost::intrusive::base_hook< z2_base_hook_matrix_set_row >
				>;

template<class Master_matrix>
class Z2_reduced_cell_set_column_with_row
		: public Z2_reduced_cell_column_with_row<Z2_set_cell,
												 Z2_set_column_type,
												 Z2_set_row_type,
												 z2_base_hook_matrix_set_row,
												 typename Master_matrix::Column_pairing_option>
{
public:
	using Cell = Z2_set_cell;
	using Column_type = Z2_set_column_type;
	using Row_type = Z2_set_row_type;

	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	Z2_reduced_cell_set_column_with_row(matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
	template<class Chain_type>
	Z2_reduced_cell_set_column_with_row(index chainIndex, const Chain_type& chain, dimension_type dimension, matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
	Z2_reduced_cell_set_column_with_row(const Z2_reduced_cell_set_column_with_row& other);

	bool is_non_zero(index rowIndex) const;

	Z2_reduced_cell_set_column_with_row& operator+=(Z2_reduced_cell_set_column_with_row& column);
	friend Z2_reduced_cell_set_column_with_row operator+(
			Z2_reduced_cell_set_column_with_row column1,
			Z2_reduced_cell_set_column_with_row const& column2){
		column1 += column2;
		return column1;
	}

private:
	using RCC = Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,z2_base_hook_matrix_set_row,typename Master_matrix::Column_pairing_option>;

	matrix_type* matrix_;
	dictionnary_type* pivotToColumnIndex_;

	void _swap_independent_rows(index rowIndex);
};

template<class Master_matrix>
inline Z2_reduced_cell_set_column_with_row<Master_matrix>::Z2_reduced_cell_set_column_with_row(matrix_type& matrix, dictionnary_type& pivotToColumnIndex)
	: RCC(), matrix_(&matrix), pivotToColumnIndex_(&pivotToColumnIndex)
{}

template<class Master_matrix>
template<class Chain_type>
inline Z2_reduced_cell_set_column_with_row<Master_matrix>::Z2_reduced_cell_set_column_with_row(
		index chainIndex, const Chain_type& chain, dimension_type dimension, matrix_type& matrix, dictionnary_type& pivotToColumnIndex)
	: RCC(chainIndex, chain, dimension), matrix_(&matrix), pivotToColumnIndex_(&pivotToColumnIndex)
{}

template<class Master_matrix>
inline Z2_reduced_cell_set_column_with_row<Master_matrix>::Z2_reduced_cell_set_column_with_row(
		const Z2_reduced_cell_set_column_with_row& other)
	: RCC(other), matrix_(other.matrix_), pivotToColumnIndex_(other.pivotToColumnIndex_)
{}

template<class Master_matrix>
inline bool Z2_reduced_cell_set_column_with_row<Master_matrix>::is_non_zero(index rowIndex) const
{
	return RCC::get_column().find(Cell(pivotToColumnIndex_->at(RCC::get_pivot()), rowIndex)) != RCC::get_column().end();
}

template<class Master_matrix>
inline Z2_reduced_cell_set_column_with_row<Master_matrix> &Z2_reduced_cell_set_column_with_row<Master_matrix>::operator+=(Z2_reduced_cell_set_column_with_row &column)
{
	Column_type& tc = RCC::get_column();
	Column_type& sc = column.get_column();
	index pos = pivotToColumnIndex_->at(RCC::get_pivot());

	for (Cell &cell : sc) {
		auto it1 = tc.find(cell);
		if (it1 != tc.end()) {
			Cell * tmp_ptr = &(*it1);
			it1->z2_base_hook_matrix_set_row::unlink();
			tc.erase(it1);
			delete tmp_ptr;
		} else {
			Cell *new_cell = new Cell(pos, cell.get_row_index());
			tc.insert(tc.end(), *new_cell);
			matrix_->at(pivotToColumnIndex_->at(cell.get_row_index())).get_row().push_back(*new_cell);//row link,no order
		}
	}

	if (!is_non_zero(RCC::get_lowest_simplex_index())){
		RCC::swap_lowest_simplex_index(column);
		_swap_independent_rows(column.get_pivot());
	}

	return *this;
}

template<class Master_matrix>
inline void Z2_reduced_cell_set_column_with_row<Master_matrix>::_swap_independent_rows(index rowIndex)
{
	std::swap(pivotToColumnIndex_->at(RCC::get_pivot()),
			  pivotToColumnIndex_->at(rowIndex));
	matrix_->at(pivotToColumnIndex_->at(RCC::get_pivot())).swap_rows(matrix_->at(pivotToColumnIndex_->at(rowIndex)));
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_CELL_LIST_COLUMN_H
