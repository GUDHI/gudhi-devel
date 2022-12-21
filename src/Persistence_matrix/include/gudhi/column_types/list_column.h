/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef LIST_COLUMN_H
#define LIST_COLUMN_H

#include <iostream>
#include <list>
#include <vector>

#include "../utilities/utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
class List_column : public Column_pairing_option, public Row_access_option
{
public:
//	using Cell = Base_cell<Field_element_type>;
	using Cell = Cell_type;
	using Column_type = std::list<Cell>;
	using iterator = typename Column_type::iterator;
	using const_iterator = typename Column_type::const_iterator;

	List_column();
	template<class Container_type>
	List_column(const Container_type& nonZeroRowIndices);
	template<class Container_type>
	List_column(const Container_type& nonZeroRowIndices, dimension_type dimension);
	template<class Row_container_type>
	List_column(index columnIndex, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	List_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	List_column(index columnIndex, const Container_type& nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer);
	List_column(const List_column& column);
	List_column(const List_column& column, index columnIndex);
	List_column(List_column&& column) noexcept;
	~List_column();

	std::vector<Field_element_type> get_content(unsigned int columnLength) const;
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	dimension_type get_dimension() const;

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;

	List_column& operator+=(List_column const &column);
	friend List_column operator+(List_column column1, List_column const& column2){
		column1 += column2;
		return column1;
	}

	List_column& operator*=(unsigned int v);
	friend List_column operator*(List_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend List_column operator*(unsigned int const& v, List_column column){
		column *= v;
		return column;
	}

	List_column& operator=(List_column other);

	friend void swap(List_column& col1, List_column& col2){
		swap(static_cast<Column_pairing_option&>(col1),
			 static_cast<Column_pairing_option&>(col2));
		swap(static_cast<Row_access_option&>(col1),
			 static_cast<Row_access_option&>(col2));
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
	}

protected:
	dimension_type dim_;
	Column_type column_;

	void _delete_cell(iterator& it);
	void _insert_cell(const Field_element_type& value, index rowIndex, const iterator& position);
	void _update_cell(const Field_element_type& value, index rowIndex, const iterator& position);
};

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_column() : dim_(0)
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_column(const Container_type &nonZeroRowIndices)
	: dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1), column_(nonZeroRowIndices.size())
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	auto it = column_.begin();
	for (const std::pair<index,Field_element_type>& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, it++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_column(const Container_type &nonZeroRowIndices, dimension_type dimension)
	: dim_(dimension), column_(nonZeroRowIndices.size())
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	auto it = column_.begin();
	for (const std::pair<index,Field_element_type>& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, it++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Row_container_type>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_column(
		index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(0)
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type, class Row_container_type>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_column(
		index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1), column_(nonZeroRowIndices.size())
{
	auto it = column_.begin();
	for (const std::pair<index,Field_element_type>& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, it++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type, class Row_container_type>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_column(
		index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(dimension), column_(nonZeroRowIndices.size())
{
	auto it = column_.begin();
	for (const std::pair<index,Field_element_type>& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, it++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_column(const List_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_)
{
	static_assert(!Row_access_option::isActive_,
			"Copy constructor not available when row access option enabled.");
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_column(
		const List_column &column, index columnIndex)
	: Column_pairing_option(column),
	  Row_access_option(columnIndex, *column.rows_),
	  dim_(column.dim_),
	  column_(column.column_.size())
{
	auto it = column_.begin();
	for (const Cell& cell : column.column_){
		_update_cell(cell.get_element(), cell.get_row_index(), it++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_column(List_column &&column) noexcept
	: Column_pairing_option(std::move(column)),
	  Row_access_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::~List_column()
{
	if constexpr (Row_access_option::isActive_){
		for (Cell& cell : column_)
			Row_access_option::unlink(&cell);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline std::vector<Field_element_type> List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::get_content(unsigned int columnLength) const
{
	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = it->get_element();
	}
	return container;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline bool List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::is_non_zero(index rowIndex) const
{
	for (const Cell& v : column_){
		if (v.get_row_index() == rowIndex) return true;
	}
	return false;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline bool List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::is_empty() const
{
	return column_.empty();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline dimension_type List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::iterator
List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::begin() noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::const_iterator
List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::iterator
List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::end() noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::const_iterator
List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::end() const noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option> &List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::operator+=(List_column const &column)
{
	if (column.is_empty()) return *this;
	if (column_.empty()){
		if constexpr (Row_access_option::isActive_){
			column_.resize(column.column_.size());
			auto it = column_.begin();
			for (const Cell& cell : column.column_)
				_update_cell(cell.get_element(), cell.get_row_index(), it++);
		} else {
			std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
		}
		return *this;
	}

	const_iterator itToAdd = column.column_.begin();
	iterator itTarget = column_.begin();
	index curRowToAdd = itToAdd->get_row_index();
	index curRowTarget = itTarget->get_row_index();

	while (itToAdd != column.column_.end() && itTarget != column_.end())
	{
		if (curRowToAdd == curRowTarget){
			itTarget->get_element() += itToAdd->get_element();
			if (itTarget->get_element() == 0u) _delete_cell(itTarget);
			else {
				if constexpr (Row_access_option::isActive_)
						Row_access_option::update_cell(*itTarget);
				itTarget++;
			}
			itToAdd++;
		} else if (curRowToAdd < curRowTarget){
			_insert_cell(itToAdd->get_element(), curRowToAdd, itTarget);
			itToAdd++;
		} else {
			itTarget++;
		}

		if (itToAdd != column.column_.end()) curRowToAdd = itToAdd->get_row_index();
		if (itTarget != column_.end()) curRowTarget = itTarget->get_row_index();
	}

	while (itToAdd != column.column_.end()){
		curRowToAdd = itToAdd->get_row_index();
		if (itToAdd != column.column_.end()){
			_insert_cell(itToAdd->get_element(), curRowToAdd, column_.end());
			itToAdd++;
		}
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option> &List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::operator*=(unsigned int v)
{
	v %= Field_element_type::get_characteristic();

	if (v == 0) {
		if constexpr (Row_access_option::isActive_){
			for (Cell& cell : column_)
				Row_access_option::unlink(&cell);
		}
		column_.clear();
		return *this;
	}

	for (Cell& cell : column_){
		cell.get_element() *= v;
		if constexpr (Row_access_option::isActive_)
				Row_access_option::update_cell(cell);
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option> &List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::operator=(List_column other)
{
	static_assert (!Row_access_option::isActive_, "= assignement not enabled with row access option.");

	Column_pairing_option::operator=(other);
	std::swap(dim_, other.dim_);
	column_.swap(other.column_);
	return *this;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline void List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::_delete_cell(iterator &it)
{
	if constexpr (Row_access_option::isActive_)
		Row_access_option::unlink(&(*it));
	column_.erase(it++);
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline void List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::_insert_cell(
		const Field_element_type &value, index rowIndex, const iterator &position)
{
	if constexpr (Row_access_option::isActive_){
		auto it = column_.insert(position, Cell(value, Row_access_option::columnIndex_, rowIndex));
		Row_access_option::insert_cell(rowIndex, &(*it));
	} else {
		column_.emplace(position, value, rowIndex);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline void List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::_update_cell(
		const Field_element_type &value, index rowIndex, const iterator &position)
{
	if constexpr (Row_access_option::isActive_){
		*position = Cell(value, Row_access_option::columnIndex_, rowIndex);
		Row_access_option::insert_cell(rowIndex, &(*position));
	} else {
		*position = Cell(value, rowIndex);
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // LIST_COLUMN_H
