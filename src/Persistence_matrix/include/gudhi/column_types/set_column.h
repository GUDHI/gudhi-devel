/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SETCOLUMN_H
#define SETCOLUMN_H

#include <iostream>
#include <list>
#include <set>

#include "../utilities/utilities.h"
#include "../utilities/Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Column_pairing_option>
class Set_column : public Column_pairing_option
{
public:
	using Cell = Base_cell<Field_element_type>;
	using iterator = typename std::set<Cell>::iterator;
	using const_iterator = typename std::set<Cell>::const_iterator;

	Set_column();
	template<class Boundary_type>
	Set_column(const Boundary_type& boundary);
	template<class Boundary_type>
	Set_column(const Boundary_type& boundary, dimension_type dimension);
//	Set_column(Set_column& column);
	Set_column(const Set_column& column);
	Set_column(Set_column&& column) noexcept;

	std::vector<Field_element_type> get_content(unsigned int columnLength) const;
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	dimension_type get_dimension() const;
	int get_pivot() const;
	Field_element_type get_pivot_value() const;
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;

	Set_column& operator+=(Set_column const &column);
	friend Set_column operator+(Set_column column1, Set_column const& column2){
		column1 += column2;
		return column1;
	}
	Set_column& operator*=(unsigned int v);
	friend Set_column operator*(Set_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Set_column operator*(unsigned int const& v, Set_column const column){
		column *= v;
		return column;
	}

	Set_column& operator=(Set_column other);

	friend void swap(Set_column& col1, Set_column& col2){
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
	}

private:
	int dim_;
	std::set<Cell> column_;
};

template<class Field_element_type, class Column_pairing_option>
inline Set_column<Field_element_type,Column_pairing_option>::Set_column() : dim_(0)
{}

template<class Field_element_type, class Column_pairing_option>
template<class Boundary_type>
inline Set_column<Field_element_type,Column_pairing_option>::Set_column(const Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1)
{
	for (const std::pair<index,Field_element_type>& p : boundary){
		column_.emplace(p.second, p.first);
	}
}

template<class Field_element_type, class Column_pairing_option>
template<class Boundary_type>
inline Set_column<Field_element_type,Column_pairing_option>::Set_column(const Boundary_type &boundary, dimension_type dimension)
	: dim_(dimension)
{
	for (const std::pair<index,Field_element_type>& p : boundary){
		column_.emplace(p.second, p.first);
	}
}

//template<class Field_element_type, class Column_pairing_option>
//inline Set_column<Field_element_type,Column_pairing_option>::Set_column(Set_column &column)
//	: Column_pairing_option(column),
//	  dim_(column.dim_),
//	  column_(column.column_)
//{}

template<class Field_element_type, class Column_pairing_option>
inline Set_column<Field_element_type,Column_pairing_option>::Set_column(const Set_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_)
{}

template<class Field_element_type, class Column_pairing_option>
inline Set_column<Field_element_type,Column_pairing_option>::Set_column(Set_column &&column) noexcept
	: Column_pairing_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Field_element_type, class Column_pairing_option>
inline std::vector<Field_element_type> Set_column<Field_element_type,Column_pairing_option>::get_content(unsigned int columnLength) const
{
	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = it->get_element();
	}
	return container;
}

template<class Field_element_type, class Column_pairing_option>
inline bool Set_column<Field_element_type,Column_pairing_option>::is_non_zero(index rowIndex) const
{
	return column_.find(Cell(0, rowIndex)) != column_.end();
}

template<class Field_element_type, class Column_pairing_option>
inline bool Set_column<Field_element_type,Column_pairing_option>::is_empty() const
{
	return column_.empty();
}

template<class Field_element_type, class Column_pairing_option>
inline dimension_type Set_column<Field_element_type,Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type, class Column_pairing_option>
inline int Set_column<Field_element_type,Column_pairing_option>::get_pivot() const
{
	if (column_.empty()) return -1;
	return column_.rbegin()->get_row_index();
}

template<class Field_element_type, class Column_pairing_option>
inline Field_element_type Set_column<Field_element_type,Column_pairing_option>::get_pivot_value() const
{
	if (column_.empty()) return 0;
	return column_.rbegin()->get_element();
}

template<class Field_element_type, class Column_pairing_option>
inline void Set_column<Field_element_type,Column_pairing_option>::clear()
{
	column_.clear();
}

template<class Field_element_type, class Column_pairing_option>
inline void Set_column<Field_element_type,Column_pairing_option>::clear(index rowIndex)
{
	column_.erase(Cell(0, rowIndex));
}

template<class Field_element_type, class Column_pairing_option>
inline void Set_column<Field_element_type,Column_pairing_option>::reorder(std::vector<index> &valueMap)
{
	std::set<Cell> newSet;
	for (const Cell& v : column_) newSet.insert(Cell(v.get_element(), valueMap[v.get_row_index()]));
	column_.swap(newSet);
}

template<class Field_element_type, class Column_pairing_option>
inline typename Set_column<Field_element_type,Column_pairing_option>::iterator
Set_column<Field_element_type,Column_pairing_option>::begin() noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Column_pairing_option>
inline typename Set_column<Field_element_type,Column_pairing_option>::const_iterator
Set_column<Field_element_type,Column_pairing_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Column_pairing_option>
inline typename Set_column<Field_element_type,Column_pairing_option>::iterator
Set_column<Field_element_type,Column_pairing_option>::end() noexcept
{
	return column_.end();
}

template<class Field_element_type, class Column_pairing_option>
inline typename Set_column<Field_element_type,Column_pairing_option>::const_iterator
Set_column<Field_element_type,Column_pairing_option>::end() const noexcept
{
	return column_.end();
}

template<class Field_element_type, class Column_pairing_option>
inline Set_column<Field_element_type,Column_pairing_option> &Set_column<Field_element_type,Column_pairing_option>::operator+=(Set_column const &column)
{
	for (const Cell& v : column.column_){
		auto c = column_.find(v);
		if (c != column_.end()){
			Cell newCell(*c);
			newCell.get_element() += v.get_element();
			column_.erase(c);
			if (newCell.get_element() != 0u) column_.insert(newCell);
		} else
			column_.insert(v);
	}

	return *this;
}

template<class Field_element_type, class Column_pairing_option>
inline Set_column<Field_element_type,Column_pairing_option> &Set_column<Field_element_type,Column_pairing_option>::operator*=(unsigned int v)
{
	v %= Field_element_type::get_characteristic();

	if (v == 0) {
		column_.clear();
		return *this;
	}

	std::set<Cell> newColumn;

	for (const Cell& cell : column_){
		Cell newCell(cell);
		newCell.get_element() *= v;
		newColumn.insert(newCell);
	}

	column_.swap(newColumn);

	return *this;
}

template<class Field_element_type, class Column_pairing_option>
inline Set_column<Field_element_type,Column_pairing_option> &Set_column<Field_element_type,Column_pairing_option>::operator=(Set_column other)
{
	std::swap(dim_, other.dim_);
	column_.swap(other.column_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // SETCOLUMN_H
