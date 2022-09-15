/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_LISTCOLUMN_H
#define Z2_LISTCOLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>

#include "../utilities/utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Column_pairing_option>
class Z2_list_column : public Column_pairing_option
{
public:
	using Cell = Z2_base_cell;
	using iterator = typename std::list<Cell>::iterator;
	using const_iterator = typename std::list<Cell>::const_iterator;

	Z2_list_column();
	template<class Boundary_type>
	Z2_list_column(Boundary_type& boundary);
	template<class Boundary_type>
	Z2_list_column(Boundary_type& boundary, dimension_type dimension);
	Z2_list_column(Z2_list_column& column);
	Z2_list_column(const Z2_list_column& column);
	Z2_list_column(Z2_list_column&& column) noexcept;

	std::vector<bool> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	dimension_type get_dimension() const;
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;

	Z2_list_column& operator+=(Z2_list_column const &column);
	template<class Friend_column_pairing_option>
	friend Z2_list_column<Friend_column_pairing_option> operator+(
			Z2_list_column<Friend_column_pairing_option> column1,
			Z2_list_column<Friend_column_pairing_option> const& column2);

	Z2_list_column& operator=(Z2_list_column other);

	template<class Friend_column_pairing_option>
	friend void swap(Z2_list_column<Friend_column_pairing_option>& col1,
					 Z2_list_column<Friend_column_pairing_option>& col2);

private:
	int dim_;
	std::list<Cell> column_;
};

template<class Column_pairing_option>
inline Z2_list_column<Column_pairing_option>::Z2_list_column() : dim_(0)
{}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_list_column<Column_pairing_option>::Z2_list_column(Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.begin(), boundary.end())
{}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_list_column<Column_pairing_option>::Z2_list_column(Boundary_type &boundary, dimension_type dimension)
	: dim_(dimension),
	  column_(boundary.begin(), boundary.end())
{}

template<class Column_pairing_option>
inline Z2_list_column<Column_pairing_option>::Z2_list_column(Z2_list_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_)
{}

template<class Column_pairing_option>
inline Z2_list_column<Column_pairing_option>::Z2_list_column(const Z2_list_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_)
{}

template<class Column_pairing_option>
inline Z2_list_column<Column_pairing_option>::Z2_list_column(Z2_list_column &&column) noexcept
	: Column_pairing_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Column_pairing_option>
inline std::vector<bool> Z2_list_column<Column_pairing_option>::get_content(unsigned int columnLength)
{
	std::vector<bool> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = 1;
	}
	return container;
}

template<class Column_pairing_option>
inline bool Z2_list_column<Column_pairing_option>::is_non_zero(index rowIndex) const
{
	for (const Cell& v : column_){
		if (v.get_row_index() == rowIndex) return true;
	}
	return false;
}

template<class Column_pairing_option>
inline bool Z2_list_column<Column_pairing_option>::is_empty() const
{
	return column_.empty();
}

template<class Column_pairing_option>
inline dimension_type Z2_list_column<Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Column_pairing_option>
inline int Z2_list_column<Column_pairing_option>::get_pivot()
{
	if (column_.empty()) return -1;

	return column_.back().get_row_index();
}

template<class Column_pairing_option>
inline void Z2_list_column<Column_pairing_option>::clear()
{
	column_.clear();
}

template<class Column_pairing_option>
inline void Z2_list_column<Column_pairing_option>::clear(index rowIndex)
{
	auto it = column_.begin();
	while (it != column_.end() && it->get_row_index() != rowIndex) it++;
	if (it != column_.end()) column_.erase(it);
}

template<class Column_pairing_option>
inline void Z2_list_column<Column_pairing_option>::reorder(std::vector<index> &valueMap)
{
	std::list<Cell>::iterator it = column_.begin();
	while (it != column_.end()) {
		it->setRowIndex(valueMap.at(it->get_row_index()));
		it++;
	}
	column_.sort();
}

template<class Column_pairing_option>
inline typename Z2_list_column<Column_pairing_option>::iterator
Z2_list_column<Column_pairing_option>::begin() noexcept
{
	return column_.begin();
}

template<class Column_pairing_option>
inline typename Z2_list_column<Column_pairing_option>::const_iterator
Z2_list_column<Column_pairing_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Column_pairing_option>
inline typename Z2_list_column<Column_pairing_option>::iterator
Z2_list_column<Column_pairing_option>::end() noexcept
{
	return column_.end();
}

template<class Column_pairing_option>
inline typename Z2_list_column<Column_pairing_option>::const_iterator
Z2_list_column<Column_pairing_option>::end() const noexcept
{
	return column_.end();
}

template<class Column_pairing_option>
inline Z2_list_column<Column_pairing_option> &Z2_list_column<Column_pairing_option>::operator+=(Z2_list_column const &column)
{
	if (column.is_empty()) return *this;
	if (column_.empty()){
		std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
		return *this;
	}

	std::list<Cell>::const_iterator itToAdd = column.column_.begin();
	std::list<Cell>::iterator itTarget = column_.begin();
	unsigned int valToAdd = itToAdd->get_row_index();
	unsigned int valTarget = itTarget->get_row_index();

	while (itToAdd != column.column_.end() && itTarget != column_.end())
	{
		if (itToAdd != column.column_.end() && itTarget != column_.end()){
			if (valToAdd == valTarget){
				column_.erase(itTarget++);
				itToAdd++;
			} else if (valToAdd < valTarget){
				column_.insert(itTarget, valToAdd);
				itToAdd++;
			} else {
				itTarget++;
			}
		}

		valToAdd = itToAdd->get_row_index();
		valTarget = itTarget->get_row_index();
	}

	while (itToAdd != column.column_.end()){
		valToAdd = itToAdd->get_row_index();
		if (itToAdd != column.column_.end()){
			column_.push_back(valToAdd);
			itToAdd++;
		}
	}

	return *this;
}

template<class Column_pairing_option>
inline Z2_list_column<Column_pairing_option> &Z2_list_column<Column_pairing_option>::operator=(Z2_list_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	return *this;
}

template<class Friend_column_pairing_option>
Z2_list_column<Friend_column_pairing_option> operator+(
		Z2_list_column<Friend_column_pairing_option> column1,
		Z2_list_column<Friend_column_pairing_option> const& column2)
{
	column1 += column2;
	return column1;
}

template<class Friend_column_pairing_option>
inline void swap(Z2_list_column<Friend_column_pairing_option>& col1,
				 Z2_list_column<Friend_column_pairing_option>& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_LISTCOLUMN_H
