/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_CELL_COLUMN_H
#define Z2_CELL_COLUMN_H

#include <unordered_map>
#include <iostream>

#include <boost/intrusive/list.hpp>
#include <boost/intrusive/set.hpp>

#include "../options.h"
#include "../utilities/utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
class Z2_reduced_cell_column_with_row : public Column_pairing_option
{
public:
	using iterator = typename Column_type::iterator;
	using const_iterator = typename Column_type::const_iterator;

	Z2_reduced_cell_column_with_row();
	template<class Chain_type>
	Z2_reduced_cell_column_with_row(index chainIndex, const Chain_type& chain, dimension_type dimension);
	Z2_reduced_cell_column_with_row(const Z2_reduced_cell_column_with_row& other);
	~Z2_reduced_cell_column_with_row();

	Column_type& get_column();
	const Column_type& get_column() const;
	Row_type& get_row();
	const Row_type& get_row() const;
	int get_pivot() const;
	int get_lowest_simplex_index() const;
	dimension_type get_dimension() const;
	bool is_empty() const;

	void swap_rows(Z2_reduced_cell_column_with_row& other);
	void swap_lowest_simplex_index(Z2_reduced_cell_column_with_row& other);

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;

private:
	Column_type column_;
	Row_type row_;
	int pivot_;		//simplex index associated to the chain
	int lowestSimplexIndex_;
	dimension_type dim_;
};

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::Z2_reduced_cell_column_with_row()
	: pivot_(-1), lowestSimplexIndex_(-1), dim_(-1)
{}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
template<class Chain_type>
inline Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::Z2_reduced_cell_column_with_row(
		index chainIndex, const Chain_type& chain, dimension_type dimension)
	: pivot_(*(chain.rbegin())), lowestSimplexIndex_(pivot_), dim_(dimension)
{
	for (index id : chain){
		Cell *new_cell = new Cell(chainIndex, id);
		column_.insert(column_.end(), *new_cell);
	}
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::Z2_reduced_cell_column_with_row(
		const Z2_reduced_cell_column_with_row& other)
	: Column_pairing_option(other),
	  pivot_(other.pivot_),
	  lowestSimplexIndex_(other.lowestSimplexIndex_),
	  dim_(other.dim_)
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

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::~Z2_reduced_cell_column_with_row()
{
	for (typename Column_type::iterator c_it = column_.begin(); c_it != column_.end(); )
	{
		auto tmp_it = c_it;
		++c_it;
		Cell *tmp_cell = &(*tmp_it);
		tmp_it->Row_base_hook::unlink();
		column_.erase(tmp_it);
		delete tmp_cell;
	}
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline Column_type& Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::get_column()
{
	return column_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline const Column_type& Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::get_column() const
{
	return column_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline Row_type& Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::get_row()
{
	return row_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline const Row_type& Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::get_row() const
{
	return row_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline int Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::get_pivot() const
{
	return pivot_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline int Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::get_lowest_simplex_index() const
{
	return lowestSimplexIndex_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline dimension_type Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook, Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline bool Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook, Column_pairing_option>::is_empty() const
{
	return column_.empty();
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline void Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::swap_rows(Z2_reduced_cell_column_with_row& other)
{
	row_.swap(other.row_);
	std::swap(pivot_, other.pivot_);
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline void Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::swap_lowest_simplex_index(Z2_reduced_cell_column_with_row& other)
{
	std::swap(lowestSimplexIndex_, other.lowestSimplexIndex_);
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline typename Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::iterator
Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::begin() noexcept
{
	return column_.begin();
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline typename Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::const_iterator
Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline typename Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::iterator
Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::end() noexcept
{
	return column_.end();
}

template<class Cell, class Column_type, class Row_type, class Row_base_hook, class Column_pairing_option>
inline typename Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::const_iterator
Z2_reduced_cell_column_with_row<Cell,Column_type,Row_type,Row_base_hook,Column_pairing_option>::end() const noexcept
{
	return column_.end();
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_CELL_COLUMN_H
