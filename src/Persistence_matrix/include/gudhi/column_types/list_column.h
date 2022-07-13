/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef LISTCOLUMN_H
#define LISTCOLUMN_H

#include <iostream>
#include <list>
#include <vector>

#include "../utilities.h"
//#include "../Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Column_pairing_option>
class List_column : public Column_pairing_option
{
public:
	using Cell = Base_cell<Field_element_type>;

	List_column();
	template<class Boundary_type>
	List_column(Boundary_type& boundary);
	template<class Boundary_type>
	List_column(Boundary_type& boundary, dimension_type dimension);
	List_column(const List_column& column);
	List_column(List_column&& column) noexcept;

	std::vector<Field_element_type> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	dimension_type get_dimension() const;
	int get_pivot();
	Field_element_type get_pivot_value();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);

	List_column& operator+=(List_column const &column);
	template<class Friend_field_element_type, class Friend_column_pairing_option>
	friend List_column<Friend_field_element_type,Friend_column_pairing_option> operator+(
			List_column<Friend_field_element_type,Friend_column_pairing_option> column1,
			List_column<Friend_field_element_type,Friend_column_pairing_option> const& column2);
	List_column& operator*=(unsigned int v);
	template<class Friend_field_element_type, class Friend_column_pairing_option>
	friend List_column<Friend_field_element_type,Friend_column_pairing_option> operator*(
			List_column<Friend_field_element_type,Friend_column_pairing_option> column,
			unsigned int const& v);
	template<class Friend_field_element_type, class Friend_column_pairing_option>
	friend List_column<Friend_field_element_type,Friend_column_pairing_option> operator*(
			unsigned int const& v,
			List_column<Friend_field_element_type,Friend_column_pairing_option> const column);

	List_column& operator=(List_column other);

	template<class Friend_field_element_type, class Friend_column_pairing_option>
	friend void swap(List_column<Friend_field_element_type,Friend_column_pairing_option>& col1,
					 List_column<Friend_field_element_type,Friend_column_pairing_option>& col2);

private:
	dimension_type dim_;
	std::list<Cell> column_;
};

template<class Field_element_type, class Column_pairing_option>
inline List_column<Field_element_type,Column_pairing_option>::List_column() : dim_(0)
{}

template<class Field_element_type, class Column_pairing_option>
template<class Boundary_type>
inline List_column<Field_element_type,Column_pairing_option>::List_column(Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1)
{
	for (std::pair<index,Field_element_type>& p : boundary){
		column_.emplace_back(p.second, p.first);
	}
}

template<class Field_element_type, class Column_pairing_option>
template<class Boundary_type>
inline List_column<Field_element_type,Column_pairing_option>::List_column(Boundary_type &boundary, dimension_type dimension)
	: dim_(dimension)
{
	for (std::pair<index,Field_element_type>& p : boundary){
		column_.emplace_back(p.second, p.first);
	}
}

template<class Field_element_type, class Column_pairing_option>
inline List_column<Field_element_type,Column_pairing_option>::List_column(const List_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_)
{}

template<class Field_element_type, class Column_pairing_option>
inline List_column<Field_element_type,Column_pairing_option>::List_column(List_column &&column) noexcept
	: Column_pairing_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Field_element_type, class Column_pairing_option>
inline std::vector<Field_element_type> List_column<Field_element_type,Column_pairing_option>::get_content(unsigned int columnLength)
{
	std::vector<Field_element_type,Column_pairing_option> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = it->element();
	}
	return container;
}

template<class Field_element_type, class Column_pairing_option>
inline bool List_column<Field_element_type,Column_pairing_option>::is_non_zero(index rowIndex) const
{
	for (Cell v : column_){
		if (v.get_row_index() == rowIndex) return true;
	}
	return false;
}

template<class Field_element_type, class Column_pairing_option>
inline bool List_column<Field_element_type,Column_pairing_option>::is_empty() const
{
	return column_.empty();
}

template<class Field_element_type, class Column_pairing_option>
inline dimension_type List_column<Field_element_type,Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type, class Column_pairing_option>
inline int List_column<Field_element_type,Column_pairing_option>::get_pivot()
{
	if (column_.empty()) return -1;

	return column_.back().get_row_index();
}

template<class Field_element_type, class Column_pairing_option>
inline Field_element_type List_column<Field_element_type,Column_pairing_option>::get_pivot_value()
{
	if (column_.empty()) return 0;

	return column_.back().get_element();
}

template<class Field_element_type, class Column_pairing_option>
inline void List_column<Field_element_type,Column_pairing_option>::clear()
{
	column_.clear();
}

template<class Field_element_type, class Column_pairing_option>
inline void List_column<Field_element_type,Column_pairing_option>::clear(index rowIndex)
{
	auto it = column_.begin();
	while (it != column_.end() && it->get_row_index() != rowIndex) it++;
	if (it != column_.end()) column_.erase(it);
}

template<class Field_element_type, class Column_pairing_option>
inline void List_column<Field_element_type,Column_pairing_option>::reorder(std::vector<index> &valueMap)
{
	typename std::list<Cell>::iterator it = column_.begin();
	while (it != column_.end()) {
		it->setRowIndex(valueMap.at(it->get_row_index()));
		it++;
	}
	column_.sort();
}

template<class Field_element_type, class Column_pairing_option>
inline List_column<Field_element_type,Column_pairing_option> &List_column<Field_element_type,Column_pairing_option>::operator+=(List_column const &column)
{
	if (column.is_empty()) return *this;
	if (column_.empty()){
		std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
		return *this;
	}

	typename std::list<Cell>::const_iterator itToAdd = column.column_.begin();
	typename std::list<Cell>::iterator itTarget = column_.begin();
	index curRowToAdd = itToAdd->get_row_index();
	index curRowTarget = itTarget->get_row_index();

	while (itToAdd != column.column_.end() && itTarget != column_.end())
	{
		if (itToAdd != column.column_.end() && itTarget != column_.end()){
			if (curRowToAdd == curRowTarget){
				itTarget->get_element() =+ itToAdd->get_element_value();
				if (itTarget->get_element() == 0u) column_.erase(itTarget++);
				else itTarget++;
				itToAdd++;
			} else if (curRowToAdd < curRowTarget){
				column_.insert(itTarget, Cell(itToAdd->get_element_value(), curRowToAdd));
				itToAdd++;
			} else {
				itTarget++;
			}
		}

		curRowToAdd = itToAdd->get_row_index();
		curRowTarget = itTarget->get_row_index();
	}

	while (itToAdd != column.column_.end()){
		curRowToAdd = itToAdd->get_row_index();
		if (itToAdd != column.column_.end()){
			column_.push_back(Cell(itToAdd->get_element_value(), curRowToAdd));
			itToAdd++;
		}
	}

	return *this;
}

template<class Field_element_type, class Column_pairing_option>
inline List_column<Field_element_type,Column_pairing_option> &List_column<Field_element_type,Column_pairing_option>::operator*=(unsigned int v)
{
	v %= Field_element_type::get_characteristic();

	if (v == 0) {
		column_.clear();
		return *this;
	}

	for (Cell& cell : column_){
		cell.get_element() *= v;
	}

	return *this;
}

template<class Field_element_type, class Column_pairing_option>
inline List_column<Field_element_type,Column_pairing_option> &List_column<Field_element_type,Column_pairing_option>::operator=(List_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	return *this;
}

template<class Friend_field_element_type, class Friend_column_pairing_option>
List_column<Friend_field_element_type,Friend_column_pairing_option> operator+(
		List_column<Friend_field_element_type,Friend_column_pairing_option> column1,
		List_column<Friend_field_element_type,Friend_column_pairing_option> const& column2)
{
	column1 += column2;
	return column1;
}

template<class Friend_field_element_type, class Friend_column_pairing_option>
List_column<Friend_field_element_type,Friend_column_pairing_option> operator*(
		List_column<Friend_field_element_type,Friend_column_pairing_option> column, unsigned int const& v)
{
	column *= v;
	return column;
}

template<class Friend_field_element_type, class Friend_column_pairing_option>
List_column<Friend_field_element_type,Friend_column_pairing_option> operator*(
		unsigned int const& v, List_column<Friend_field_element_type,Friend_column_pairing_option> column)
{
	column *= v;
	return column;
}

template<class Friend_field_element_type, class Friend_column_pairing_option>
inline void swap(List_column<Friend_field_element_type,Friend_column_pairing_option>& col1,
				 List_column<Friend_field_element_type,Friend_column_pairing_option>& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // LISTCOLUMN_H
