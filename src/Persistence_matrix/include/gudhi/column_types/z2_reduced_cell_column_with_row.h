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
#include "../utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<Column_types boost_column_type>
class Z2_reduced_cell_column_with_row
{
public:
	struct Cell;

private:
	struct matrix_row_tag;
	struct matrix_column_tag;

	using base_hook_list = boost::intrusive::list_base_hook <
				boost::intrusive::tag < matrix_column_tag >
			  , boost::intrusive::link_mode < boost::intrusive::safe_link >
		>;
	using base_hook_set = boost::intrusive::set_base_hook <
				boost::intrusive::tag < matrix_column_tag >
			  , boost::intrusive::optimize_size<true>
			  , boost::intrusive::link_mode < boost::intrusive::safe_link >
		>;

	using column_list = boost::intrusive::list <
				Cell
			  , boost::intrusive::constant_time_size<false>
			  , boost::intrusive::base_hook< base_hook_list >  >;
	using column_set = boost::intrusive::set <
				Cell
			  , boost::intrusive::constant_time_size<false>
			  , boost::intrusive::base_hook< base_hook_set >  >;

	using base_hook_matrix_row = boost::intrusive::list_base_hook<
				boost::intrusive::tag < matrix_row_tag >
			  , boost::intrusive::link_mode < boost::intrusive::auto_unlink >
		>;
	using base_hook_matrix_column = typename std::conditional<
				boost_column_type == Column_types::SET,
				base_hook_set,
				base_hook_list >::type;

public:
	struct Cell : public Z2_row_cell, public base_hook_matrix_row, public base_hook_matrix_column
	{
		Cell(index columnIndex, index rowIndex)
			: Z2_row_cell(columnIndex, rowIndex){};
	};

	using Column_type = typename std::conditional<
						boost_column_type == Column_types::SET,
						column_set,
						column_list
					>::type;
	using Row_type = boost::intrusive::list <
						Cell
					  , boost::intrusive::constant_time_size<false>
					  , boost::intrusive::base_hook< base_hook_matrix_row >
					>;

	Z2_reduced_cell_column_with_row();
	template<class Chain_type>
	Z2_reduced_cell_column_with_row(index chainIndex, Chain_type& chain);
	template<class Chain_type>
	Z2_reduced_cell_column_with_row(index chainIndex, Chain_type& chain, index pairedColumnIndex);
	Z2_reduced_cell_column_with_row(const Z2_reduced_cell_column_with_row& other);
	~Z2_reduced_cell_column_with_row();

	index get_paired_chain_index();
	bool is_paired();
	void assign_paired_chain(index other_col);
	void unassign_paired_chain();
	Column_type& get_column();
	Row_type& get_row();
	int get_pivot();
	int get_lowest_simplex_index();
	void swap_rows(Z2_reduced_cell_column_with_row& other);
	void swap_lowest_simplex_index(Z2_reduced_cell_column_with_row& other);

private:
	Column_type column_;
	Row_type row_;
	int pivot_;		//simplex index associated to the chain
	int lowestSimplexIndex_;
	int pairedColumn_;
};

template<Column_types boost_column_type>
inline Z2_reduced_cell_column_with_row<boost_column_type>::Z2_reduced_cell_column_with_row()
	: pivot_(-1), lowestSimplexIndex_(-1), pairedColumn_(-1)
{}

template<Column_types boost_column_type>
template<class Chain_type>
inline Z2_reduced_cell_column_with_row<boost_column_type>::Z2_reduced_cell_column_with_row(
		index chainIndex, Chain_type& chain)
	: pivot_(*(chain.rbegin())), lowestSimplexIndex_(pivot_), pairedColumn_(-1)
{
	for (index id : chain){
		Cell *new_cell = new Cell(chainIndex, id);
		column_.insert(column_.end(), *new_cell);
	}
}

template<Column_types boost_column_type>
template<class Chain_type>
inline Z2_reduced_cell_column_with_row<boost_column_type>::Z2_reduced_cell_column_with_row(
		index chainIndex, Chain_type& chain, index pairedColumnIndex)
	: pivot_(*(chain.rbegin())), lowestSimplexIndex_(pivot_), pairedColumn_(pairedColumnIndex)
{
	for (index id : chain){
		Cell *new_cell = new Cell(chainIndex, id);
		column_.insert(column_.end(), *new_cell);
	}
}

template<Column_types boost_column_type>
inline Z2_reduced_cell_column_with_row<boost_column_type>::Z2_reduced_cell_column_with_row(
		const Z2_reduced_cell_column_with_row& other)
	: pivot_(other.pivot_),
	  lowestSimplexIndex_(other.lowestSimplexIndex_),
	  pairedColumn_(other.pairedColumn_)
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

template<Column_types boost_column_type>
inline Z2_reduced_cell_column_with_row<boost_column_type>::~Z2_reduced_cell_column_with_row()
{ //empty the column, call delete on all cells
	for (typename Column_type::iterator c_it = column_.begin(); c_it != column_.end(); )
	{
		auto tmp_it = c_it;
		++c_it;
		Cell *tmp_cell = &(*tmp_it);
		tmp_it->base_hook_matrix_row::unlink(); //rm from row
		column_.erase(tmp_it);
		delete tmp_cell;
	}
}

template<Column_types boost_column_type>
inline index Z2_reduced_cell_column_with_row<boost_column_type>::get_paired_chain_index()
{
	assert(pairedColumn_ != -1 && "Column not paired.");
	return static_cast<index>(pairedColumn_);
}

template<Column_types boost_column_type>
inline bool Z2_reduced_cell_column_with_row<boost_column_type>::is_paired()
{
	return pairedColumn_ != -1;
}

template<Column_types boost_column_type>
inline void Z2_reduced_cell_column_with_row<boost_column_type>::assign_paired_chain(
		index other_col)
{
	pairedColumn_ = other_col;
}

template<Column_types boost_column_type>
inline void Z2_reduced_cell_column_with_row<boost_column_type>::unassign_paired_chain()
{
	pairedColumn_ = -1;
}

template<Column_types boost_column_type>
inline typename Z2_reduced_cell_column_with_row<boost_column_type>::Column_type&
Z2_reduced_cell_column_with_row<boost_column_type>::get_column()
{
	return column_;
}

template<Column_types boost_column_type>
inline typename Z2_reduced_cell_column_with_row<boost_column_type>::Row_type&
Z2_reduced_cell_column_with_row<boost_column_type>::get_row()
{
	return row_;
}

template<Column_types boost_column_type>
inline int Z2_reduced_cell_column_with_row<boost_column_type>::get_pivot()
{
	return pivot_;
}

template<Column_types boost_column_type>
inline int Z2_reduced_cell_column_with_row<boost_column_type>::get_lowest_simplex_index()
{
	return lowestSimplexIndex_;
}

template<Column_types boost_column_type>
inline void Z2_reduced_cell_column_with_row<boost_column_type>::swap_rows(Z2_reduced_cell_column_with_row& other)
{
	std::swap(row_, other.row_);
	std::swap(pivot_, other.pivot_);
}

template<Column_types boost_column_type>
inline void Z2_reduced_cell_column_with_row<boost_column_type>::swap_lowest_simplex_index(Z2_reduced_cell_column_with_row& other)
{
	std::swap(lowestSimplexIndex_, other.lowestSimplexIndex_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_CELL_COLUMN_H
