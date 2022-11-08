/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CELL_SET_COLUMN_H
#define CELL_SET_COLUMN_H

#include "reduced_cell_column_with_row.h"

namespace Gudhi {
namespace persistence_matrix {

struct matrix_set_row_tag;
struct matrix_set_column_tag;

using base_hook_matrix_set_row = boost::intrusive::list_base_hook<
			boost::intrusive::tag < matrix_set_row_tag >
		  , boost::intrusive::link_mode < boost::intrusive::auto_unlink >
	>;
using base_hook_matrix_set_column = boost::intrusive::set_base_hook <
			boost::intrusive::tag < matrix_set_column_tag >
		  , boost::intrusive::link_mode < boost::intrusive::safe_link >
		>;

template<class Field_element_type>
struct Set_cell : public Row_cell<Field_element_type>, public base_hook_matrix_set_row, public base_hook_matrix_set_column
{
	Set_cell(unsigned int element, index columnIndex, index rowIndex)
		: Row_cell<Field_element_type>(element, columnIndex, rowIndex){};
};

template<class Field_element_type>
using Set_column_type = boost::intrusive::set <
					Set_cell<Field_element_type>
				  , boost::intrusive::constant_time_size<false>
				  , boost::intrusive::base_hook< base_hook_matrix_set_column >  >;

template<class Field_element_type>
using Set_row_type = boost::intrusive::list <
					Set_cell<Field_element_type>
				  , boost::intrusive::constant_time_size<false>
				  , boost::intrusive::base_hook< base_hook_matrix_set_row >
				>;

template<class Master_matrix, class Field_element_type, class Column_pairing_option>
class Reduced_cell_set_column_with_row
		: public Reduced_cell_column_with_row<Set_cell<Field_element_type>,
											  Set_column_type<Field_element_type>,
											  Set_row_type<Field_element_type>,
											  base_hook_matrix_set_row,
											  Field_element_type,
											  Column_pairing_option>
{
public:
	using Cell = Set_cell<Field_element_type>;
	using Column_type = Set_column_type<Field_element_type>;
	using Row_type = Set_row_type<Field_element_type>;

	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	Reduced_cell_set_column_with_row(matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
	template<class Chain_type>
	Reduced_cell_set_column_with_row(index chainIndex, const Chain_type& chain, dimension_type dimension, matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
	Reduced_cell_set_column_with_row(const Reduced_cell_set_column_with_row& other);

	bool is_non_zero(index rowIndex) const;

	Reduced_cell_set_column_with_row& operator+=(Reduced_cell_set_column_with_row &column);
	friend Reduced_cell_set_column_with_row operator+(
			Reduced_cell_set_column_with_row column1,
			Reduced_cell_set_column_with_row const& column2){
		column1 += column2;
		return column1;
	}
	Reduced_cell_set_column_with_row& operator*=(unsigned int v);
	friend Reduced_cell_set_column_with_row operator*(
			Reduced_cell_set_column_with_row column,
			unsigned int const& v){
		column *= v;
		return column;
	}
	friend Reduced_cell_set_column_with_row operator*(
			unsigned int const& v,
			Reduced_cell_set_column_with_row column){
		column *= v;
		return column;
	}

private:
	using RCC = Reduced_cell_column_with_row<Cell,Column_type,Row_type,base_hook_matrix_set_row,typename Master_matrix::Field_type,typename Master_matrix::Column_pairing_option>;

	matrix_type* matrix_;
	dictionnary_type* pivotToColumnIndex_;

	void _swap_independent_rows(index rowIndex);
};

template<class Master_matrix, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_set_column_with_row<Master_matrix,Field_element_type,Column_pairing_option>::Reduced_cell_set_column_with_row(matrix_type& matrix, dictionnary_type& pivotToColumnIndex)
	: RCC(), matrix_(&matrix), pivotToColumnIndex_(&pivotToColumnIndex)
{}

template<class Master_matrix, class Field_element_type, class Column_pairing_option>
template<class Chain_type>
inline Reduced_cell_set_column_with_row<Master_matrix,Field_element_type,Column_pairing_option>::Reduced_cell_set_column_with_row(
		index chainIndex, const Chain_type& chain, dimension_type dimension, matrix_type& matrix, dictionnary_type& pivotToColumnIndex)
	: RCC(chainIndex, chain, dimension), matrix_(&matrix), pivotToColumnIndex_(&pivotToColumnIndex)
{}

template<class Master_matrix, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_set_column_with_row<Master_matrix,Field_element_type,Column_pairing_option>::Reduced_cell_set_column_with_row(
		const Reduced_cell_set_column_with_row& other)
	: RCC(other), matrix_(other.matrix_), pivotToColumnIndex_(other.pivotToColumnIndex_)
{}

template<class Master_matrix, class Field_element_type, class Column_pairing_option>
inline void Reduced_cell_set_column_with_row<Master_matrix,Field_element_type,Column_pairing_option>::_swap_independent_rows(index rowIndex)
{
	std::swap(pivotToColumnIndex_->at(RCC::get_pivot()),
			  pivotToColumnIndex_->at(rowIndex));
	matrix_->at(pivotToColumnIndex_->at(RCC::get_pivot())).swap_rows(matrix_->at(pivotToColumnIndex_->at(rowIndex)));
}

template<class Master_matrix, class Field_element_type, class Column_pairing_option>
inline bool Reduced_cell_set_column_with_row<Master_matrix,Field_element_type,Column_pairing_option>::is_non_zero(index rowIndex) const
{
	return RCC::get_column().find(Cell(pivotToColumnIndex_->at(RCC::get_pivot()), rowIndex)) != RCC::get_column().end();
}

template<class Master_matrix, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_set_column_with_row<Master_matrix,Field_element_type,Column_pairing_option> &Reduced_cell_set_column_with_row<Master_matrix,Field_element_type,Column_pairing_option>::operator+=(Reduced_cell_set_column_with_row &column)
{
	Column_type& tc = RCC::get_column();
	Column_type& sc = column.get_column();
	index pos = pivotToColumnIndex_->at(RCC::get_pivot());

	for (Cell &cell : sc) {
		auto it1 = tc.find(cell);
		if (it1 != tc.end()) {
			it1->get_element() += cell.get_element();
			if (it1->get_element() == 0u){
				Cell *tmp_ptr = &(*it1);
				it1->base_hook_matrix_set_row::unlink();
				tc.erase(it1);
				delete tmp_ptr;
			}
		} else {
			Cell *new_cell = new Cell(cell.get_element(), pos, cell.get_row_index());
			tc.insert(tc.end(), *new_cell);
			matrix_->at(pivotToColumnIndex_->at(cell.get_row_index())).get_row().push_back(*new_cell);
		}
	}

	if (!is_non_zero(RCC::get_lowest_simplex_index())){
		RCC::swap_lowest_simplex_index(column);
		_swap_independent_rows(column.get_pivot());
	}

	return *this;
}

template<class Master_matrix, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_set_column_with_row<Master_matrix,Field_element_type,Column_pairing_option> &Reduced_cell_set_column_with_row<Master_matrix,Field_element_type,Column_pairing_option>::operator*=(unsigned int v)
{
	v %= Field_element_type::get_characteristic();

	if (v == 0) {
		auto it = RCC::get_column().begin();
		while (it != RCC::get_column().end()){
			typename Column_type::iterator tmp_it = it;
			++it;
			Cell* tmp_ptr = &(*tmp_it);
			tmp_it->base_hook_matrix_set_row::unlink();
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

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CELL_SET_COLUMN_H
