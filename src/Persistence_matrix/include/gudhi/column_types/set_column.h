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

#include "../utilities.h"
#include "../Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Column_pairing_option>
class Set_column : public Column_pairing_option
{
public:
	using Cell = Base_cell<Field_element_type>;

	Set_column();
	template<class Boundary_type>
	Set_column(Boundary_type& boundary);
	template<class Boundary_type>
	Set_column(Boundary_type& boundary, dimension_type dimension);
	Set_column(Set_column& column);
	Set_column(const Set_column& column);
	Set_column(Set_column&& column) noexcept;

	std::vector<Field_element_type> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	Field_element_type get_pivot_value();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);

	Set_column& operator+=(Set_column const &column);
	template<class Friend_field_element_type, class Friend_column_pairing_option>
	friend Set_column<Friend_field_element_type,Friend_column_pairing_option> operator+(
			Set_column<Friend_field_element_type,Friend_column_pairing_option> column1,
			Set_column<Friend_field_element_type,Friend_column_pairing_option> const& column2);
	Set_column& operator*=(unsigned int v);
	template<class Friend_field_element_type, class Friend_column_pairing_option>
	friend Set_column<Friend_field_element_type,Friend_column_pairing_option> operator*(
			Set_column<Friend_field_element_type,Friend_column_pairing_option> column,
			unsigned int const& v);
	template<class Friend_field_element_type, class Friend_column_pairing_option>
	friend Set_column<Friend_field_element_type,Friend_column_pairing_option> operator*(
			unsigned int const& v,
			Set_column<Friend_field_element_type,Friend_column_pairing_option> const column);

	Set_column& operator=(Set_column other);

	template<class Friend_field_element_type, class Friend_column_pairing_option>
	friend void swap(Set_column<Friend_field_element_type,Friend_column_pairing_option>& col1,
					 Set_column<Friend_field_element_type,Friend_column_pairing_option>& col2);

private:
	int dim_;
	std::set<Cell> column_;
};

template<class Field_element_type, class Column_pairing_option>
inline Set_column<Field_element_type,Column_pairing_option>::Set_column() : dim_(0)
{}

template<class Field_element_type, class Column_pairing_option>
template<class Boundary_type>
inline Set_column<Field_element_type,Column_pairing_option>::Set_column(Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1)
{
	for (std::pair<index,Field_element_type>& p : boundary){
		column_.emplace(p.second, p.first);
	}
}

template<class Field_element_type, class Column_pairing_option>
template<class Boundary_type>
inline Set_column<Field_element_type,Column_pairing_option>::Set_column(Boundary_type &boundary, dimension_type dimension)
	: dim_(dimension)
{
	for (std::pair<index,Field_element_type>& p : boundary){
		column_.emplace(p.second, p.first);
	}
}

template<class Field_element_type, class Column_pairing_option>
inline Set_column<Field_element_type,Column_pairing_option>::Set_column(Set_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_)
{}

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
inline std::vector<Field_element_type> Set_column<Field_element_type,Column_pairing_option>::get_content(unsigned int columnLength)
{
	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = it->get_element_value();
	}
	return container;
}

template<class Field_element_type, class Column_pairing_option>
inline bool Set_column<Field_element_type,Column_pairing_option>::is_non_zero(index rowIndex) const
{
	return column_.find(Cell(0, rowIndex)) != column_.end();
}

template<class Field_element_type, class Column_pairing_option>
inline bool Set_column<Field_element_type,Column_pairing_option>::is_empty()
{
	return column_.empty();
}

template<class Field_element_type, class Column_pairing_option>
inline dimension_type Set_column<Field_element_type,Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type, class Column_pairing_option>
inline int Set_column<Field_element_type,Column_pairing_option>::get_pivot()
{
	if (column_.empty()) return -1;
	return column_.rbegin()->get_row_index();
}

template<class Field_element_type, class Column_pairing_option>
inline Field_element_type Set_column<Field_element_type,Column_pairing_option>::get_pivot_value()
{
	if (column_.empty()) return 0;
	Cell c = *(column_.rbegin());
	return c.get_element();
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
	for (const Cell& v : column_) newSet.insert(Cell(v.get_element(), valueMap.at(v.get_row_index())));
	column_.swap(newSet);
}

template<class Field_element_type, class Column_pairing_option>
inline Set_column<Field_element_type,Column_pairing_option> &Set_column<Field_element_type,Column_pairing_option>::operator+=(Set_column const &column)
{
	for (const Cell& v : column.column_){
		auto c = column_.find(v);
		if (c != column_.end()){
			Cell newCell(*c);
			newCell.get_element() += v.get_element_value();
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
	std::swap(column_, other.column_);
	return *this;
}

template<class Friend_field_element_type, class Friend_column_pairing_option>
Set_column<Friend_field_element_type,Friend_column_pairing_option> operator+(
		Set_column<Friend_field_element_type,Friend_column_pairing_option> column1,
		Set_column<Friend_field_element_type,Friend_column_pairing_option> const& column2)
{
	column1 += column2;
	return column1;
}

template<class Friend_field_element_type, class Friend_column_pairing_option>
Set_column<Friend_field_element_type,Friend_column_pairing_option> operator*(
		Set_column<Friend_field_element_type,Friend_column_pairing_option> column, unsigned int const& v)
{
	column *= v;
	return column;
}

template<class Friend_field_element_type, class Friend_column_pairing_option>
Set_column<Friend_field_element_type,Friend_column_pairing_option> operator*(
		unsigned int const& v, Set_column<Friend_field_element_type,Friend_column_pairing_option> column)
{
	column *= v;
	return column;
}

template<class Friend_field_element_type, class Friend_column_pairing_option>
inline void swap(Set_column<Friend_field_element_type,Friend_column_pairing_option>& col1,
				 Set_column<Friend_field_element_type,Friend_column_pairing_option>& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // SETCOLUMN_H