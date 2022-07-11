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

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
class Reduced_cell_column_with_row : Column_pairing_option
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
	struct Cell : public Row_cell<Field_element_type>, public base_hook_matrix_row, public base_hook_matrix_column
	{
		Cell(unsigned int element, index columnIndex, index rowIndex)
			: Row_cell<Field_element_type>(element, columnIndex, rowIndex){};
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

	Reduced_cell_column_with_row();
	template<class Chain_type>
	Reduced_cell_column_with_row(index chainIndex, Chain_type& chain, dimension_type dimension);
	Reduced_cell_column_with_row(const Reduced_cell_column_with_row& other);
	~Reduced_cell_column_with_row();

	Column_type& get_column();
	Row_type& get_row();
	int get_pivot();
	int get_lowest_simplex_index();
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

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::Reduced_cell_column_with_row()
	: pivot_(-1), lowestSimplexIndex_(-1), dim_(-1)
{}

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
template<class Chain_type>
inline Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::Reduced_cell_column_with_row(
		index chainIndex, Chain_type& chain, dimension_type dimension)
	: pivot_(chain.rbegin()->first), lowestSimplexIndex_(pivot_), dim_(dimension)
{
	for (std::pair<index,Field_element_type>& p : chain){
		Cell *new_cell = new Cell(p.second, chainIndex, p.first);
		column_.insert(column_.end(), *new_cell);
	}
}

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::Reduced_cell_column_with_row(
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

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::~Reduced_cell_column_with_row()
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

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline typename Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::Column_type&
Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::get_column()
{
	return column_;
}

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline typename Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::Row_type&
Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::get_row()
{
	return row_;
}

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline int Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::get_pivot()
{
	return pivot_;
}

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline int Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::get_lowest_simplex_index()
{
	return lowestSimplexIndex_;
}

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline dimension_type Reduced_cell_column_with_row<boost_column_type, Field_element_type, Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline bool Reduced_cell_column_with_row<boost_column_type, Field_element_type, Column_pairing_option>::is_empty() const
{
	return column_.empty();
}

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline void Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::swap_rows(Reduced_cell_column_with_row& other)
{
	std::swap(row_, other.row_);
	std::swap(pivot_, other.pivot_);
}

template<Column_types boost_column_type, class Field_element_type, class Column_pairing_option>
inline void Reduced_cell_column_with_row<boost_column_type,Field_element_type,Column_pairing_option>::swap_lowest_simplex_index(Reduced_cell_column_with_row& other)
{
	std::swap(lowestSimplexIndex_, other.lowestSimplexIndex_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CELL_COLUMN_H
