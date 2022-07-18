/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CELL_COLUMN_H
#define CELL_COLUMN_H

#include <unordered_map>
#include <iostream>

#include <boost/intrusive/list.hpp>
#include <boost/intrusive/set.hpp>

#include "../options.h"
#include "../utilities.h"
#include "../Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
class Reduced_cell_column_with_row : public Column_pairing_option
{
public:
	Reduced_cell_column_with_row();
	template<class Chain_type>
	Reduced_cell_column_with_row(index chainIndex, Chain_type& chain, dimension_type dimension);
	Reduced_cell_column_with_row(const Reduced_cell_column_with_row& other);
	~Reduced_cell_column_with_row();

	Column_type& get_column();
	Row_type& get_row();
	int get_pivot() const;
	Field_element_type get_pivot_value();
	int get_lowest_simplex_index() const;
	dimension_type get_dimension() const;
	bool is_empty() const;

	void swap_rows(Reduced_cell_column_with_row& other);
	void swap_lowest_simplex_index(Reduced_cell_column_with_row& other);

private:
	Column_type column_;
	Row_type row_;
	int pivot_;		//simplex index associated to the chain
	int lowestSimplexIndex_;
	dimension_type dim_;
};

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::Reduced_cell_column_with_row()
	: pivot_(-1), lowestSimplexIndex_(-1), dim_(-1)
{}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
template<class Chain_type>
inline Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::Reduced_cell_column_with_row(
		index chainIndex, Chain_type& chain, dimension_type dimension)
	: pivot_(chain.rbegin()->first), lowestSimplexIndex_(pivot_), dim_(dimension)
{
	for (const std::pair<index,Field_element_type>& p : chain){
		Cell *new_cell = new Cell(p.second, chainIndex, p.first);
		column_.insert(column_.end(), *new_cell);
	}
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::Reduced_cell_column_with_row(
		const Reduced_cell_column_with_row& other)
	: Column_pairing_option(other),
	  pivot_(other.pivot_),
	  lowestSimplexIndex_(other.lowestSimplexIndex_)
{
	//Cloner object function
	struct new_cloner
	{
	   Cell *operator()(const Cell &clone_this)
	   {  return new Cell(clone_this);  }
	};

	//The disposer object function
	struct delete_disposer
	{
	   void operator()(Cell *delete_this)
	   {  delete delete_this;  }
	};

	column_.clone_from(other.column_, new_cloner(), delete_disposer());
	row_.clone_from(other.row_, new_cloner(), delete_disposer());
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::~Reduced_cell_column_with_row()
{ //empty the column, call delete on all cells
	for (typename Column_type::iterator c_it = column_.begin(); c_it != column_.end(); )
	{
		auto tmp_it = c_it;
		++c_it;
		Cell *tmp_cell = &(*tmp_it);
		tmp_it->Row_base_hook::unlink(); //rm from row
		column_.erase(tmp_it);
		delete tmp_cell;
	}
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline Column_type& Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::get_column()
{
	return column_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline Row_type& Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::get_row()
{
	return row_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline int Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::get_pivot() const
{
	return pivot_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline Field_element_type Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::get_pivot_value()
{
	for (Cell& cell : column_){
		if (cell.get_row_index() == pivot_) return cell.get_element();
	}
	return Field_element_type();
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline int Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::get_lowest_simplex_index() const
{
	return lowestSimplexIndex_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline dimension_type Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook, Field_element_type, Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline bool Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook, Field_element_type, Column_pairing_option>::is_empty() const
{
	return column_.empty();
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline void Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::swap_rows(Reduced_cell_column_with_row& other)
{
	std::swap(row_, other.row_);
	std::swap(pivot_, other.pivot_);
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Field_element_type, class Column_pairing_option>
inline void Reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Field_element_type,Column_pairing_option>::swap_lowest_simplex_index(Reduced_cell_column_with_row& other)
{
	std::swap(lowestSimplexIndex_, other.lowestSimplexIndex_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CELL_COLUMN_H
