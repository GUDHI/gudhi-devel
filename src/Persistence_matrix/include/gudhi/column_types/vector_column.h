/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef VECTOR_COLUMN_H
#define VECTOR_COLUMN_H

#include <iostream>
#include <vector>
#include <unordered_set>

#include <boost/iterator/indirect_iterator.hpp>

#include "../utilities/utilities.h"
#include "../utilities/Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
class Vector_column : public Column_pairing_option, public Row_access_option
{
public:
//	using Cell = Base_cell<Field_element_type>;
	using Cell = Cell_type;
	using Column_type = std::vector<Cell*>;
	using iterator = boost::indirect_iterator<typename Column_type::iterator>;
	using const_iterator = boost::indirect_iterator<typename Column_type::const_iterator>;
	using reverse_iterator = typename Column_type::reverse_iterator;
	using const_reverse_iterator = typename Column_type::const_reverse_iterator;

	Vector_column();
	template<class Container_type>
	Vector_column(const Container_type& nonZeroRowIndices);
	template<class Container_type>
	Vector_column(const Container_type& nonZeroRowIndices, dimension_type dimension);
	template<class Row_container_type>
	Vector_column(index columnIndex, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Vector_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Vector_column(index columnIndex, const Container_type& nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer);
	Vector_column(const Vector_column& column);
	Vector_column(const Vector_column& column, index columnIndex);
	Vector_column(Vector_column&& column) noexcept;
	~Vector_column();

	std::vector<Field_element_type> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	dimension_type get_dimension() const;

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;
	reverse_iterator rbegin() noexcept;
	const_reverse_iterator rbegin() const noexcept;
	reverse_iterator rend() noexcept;
	const_reverse_iterator rend() const noexcept;

	Vector_column& operator+=(Vector_column const &column);
	friend Vector_column operator+(Vector_column column1, Vector_column const &column2){
		column1 += column2;
		return column1;
	}
	Vector_column& operator*=(unsigned int v);
	friend Vector_column operator*(Vector_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Vector_column operator*(unsigned int const& v, Vector_column column){
		column *= v;
		return column;
	}

	Vector_column& operator=(Vector_column other);

	friend void swap(Vector_column& col1, Vector_column& col2){
		swap(static_cast<Column_pairing_option&>(col1),
			 static_cast<Column_pairing_option&>(col2));
		swap(static_cast<Row_access_option&>(col1),
			 static_cast<Row_access_option&>(col2));
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
	}

protected:
	using real_iterator = typename Column_type::iterator;
	using real_const_iterator = typename Column_type::const_iterator;

	int dim_;
	Column_type column_;

	void _delete_cell(Cell* cell);
	void _insert_cell(const Field_element_type& value, index rowIndex, Column_type& column);
	void _update_cell(const Field_element_type& value, index rowIndex, index position);
};

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::Vector_column() : dim_(0)
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::Vector_column(const Container_type &nonZeroRowIndices)
	: dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
	  column_(nonZeroRowIndices.size(), nullptr)
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	unsigned int i = 0;
	for (const std::pair<index,Field_element_type>& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, i++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::Vector_column(const Container_type &nonZeroRowIndices, dimension_type dimension)
	: dim_(dimension),
	  column_(nonZeroRowIndices.size(), nullptr)
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	unsigned int i = 0;
	for (const std::pair<index,Field_element_type>& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, i++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Row_container_type>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::Vector_column(
		index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(0)
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::Vector_column(
		index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1), column_(nonZeroRowIndices.size(), nullptr)
{
	unsigned int i = 0;
	for (const std::pair<index,Field_element_type>& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, i++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::Vector_column(
		index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(dimension), column_(nonZeroRowIndices.size(), nullptr)
{
	unsigned int i = 0;
	for (const std::pair<index,Field_element_type>& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, i++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::Vector_column(const Vector_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_.size(), nullptr)
{
	static_assert(!Row_access_option::isActive_,
			"Copy constructor not available when row access option enabled.");

	unsigned int i = 0;
	for (const Cell* cell : column.column_){
		_update_cell(cell->get_element(), cell->get_row_index(), i++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::Vector_column(
		const Vector_column &column, index columnIndex)
	: Column_pairing_option(column),
	  Row_access_option(columnIndex, *column.rows_),
	  dim_(column.dim_),
	  column_(column.column_.size(), nullptr)
{
	unsigned int i = 0;
	for (const Cell* cell : column.column_){
		_update_cell(cell->get_element(), cell->get_row_index(), i++);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::Vector_column(Vector_column &&column) noexcept
	: Column_pairing_option(std::move(column)),
	  Row_access_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::~Vector_column()
{
	for (Cell* cell : column_){
		_delete_cell(cell);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline std::vector<Field_element_type> Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::get_content(unsigned int columnLength)
{
	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && (*it)->get_row_index() < columnLength; ++it){
		container[(*it)->get_row_index()] = (*it)->get_element();
	}
	return container;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline bool Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::is_non_zero(index rowIndex) const
{
	for (const Cell* v : column_){
		if (v->get_row_index() == rowIndex) return true;
	}
	return false;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline bool Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::is_empty() const
{
	return column_.empty();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline dimension_type Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::iterator
Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::begin() noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::const_iterator
Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::iterator
Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::end() noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::const_iterator
Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::end() const noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::reverse_iterator
Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::rbegin() noexcept
{
	return column_.rbegin();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::const_reverse_iterator
Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::rbegin() const noexcept
{
	return column_.rbegin();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::reverse_iterator
Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::rend() noexcept
{
	return column_.rend();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::const_reverse_iterator
Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::rend() const noexcept
{
	return column_.rend();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option> &Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::operator+=(const Vector_column &column)
{
	if (column.is_empty()) return *this;

	if (column_.empty()){
		if constexpr (Row_access_option::isActive_){
			column_.resize(column.column_.size());
			unsigned int i = 0;
			for (const Cell* cell : column.column_)
				_update_cell(cell->get_element(), cell->get_row_index(), i++);
		} else {
			std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
		}
		return *this;
	}

	Column_type newColumn;

	real_const_iterator itToAdd = column.column_.begin();
	real_iterator itTarget = column_.begin();

	while (itToAdd != column.column_.end() && itTarget != column_.end())
	{
		Cell* cellToAdd = *itToAdd;
		Cell* cellTarget = *itTarget;
		unsigned int curRowToAdd = cellToAdd->get_row_index();
		unsigned int curRowTarget = cellTarget->get_row_index();

		if (curRowToAdd == curRowTarget){
			Field_element_type sum = cellTarget->get_element() + cellToAdd->get_element();
			if (sum != 0) {
				cellTarget->set_element(sum);
				newColumn.push_back(cellTarget);
				if constexpr (Row_access_option::isActive_)
						Row_access_option::update_cell(*cellTarget);
			} else {
				_delete_cell(cellTarget);
			}
			itTarget++;
			itToAdd++;
		} else if (curRowToAdd < curRowTarget){
			_insert_cell(cellToAdd->get_element(), curRowToAdd, newColumn);
			itToAdd++;
		} else {
			newColumn.push_back(cellTarget);
			itTarget++;
		}
	}

	while (itToAdd != column.column_.end()){
		_insert_cell((*itToAdd)->get_element(), (*itToAdd)->get_row_index(), newColumn);
		itToAdd++;
	}

	while (itTarget != column_.end()){
		newColumn.push_back(*itTarget);
		itTarget++;
	}

	column_.swap(newColumn);

	return *this;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option> &Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::operator*=(unsigned int v)
{
	v %= Field_element_type::get_characteristic();

	if (v == 0) {
		for (Cell* cell : column_){
			_delete_cell(cell);
		}
		column_.clear();
		return *this;
	}

	for (Cell* cell : column_){
		cell->get_element() *= v;
		if constexpr (Row_access_option::isActive_)
				Row_access_option::update_cell(*cell);
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option> &Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::operator=(Vector_column other)
{
	static_assert (!Row_access_option::isActive_, "= assignement not enabled with row access option.");

	Column_pairing_option::operator=(other);
	std::swap(dim_, other.dim_);
	column_.swap(other.column_);
	return *this;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline void Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::_delete_cell(Cell* cell)
{
	if constexpr (Row_access_option::isActive_)
		Row_access_option::unlink(cell);
	delete cell;
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline void Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::_insert_cell(const Field_element_type &value, index rowIndex, Column_type &column)
{
	if constexpr (Row_access_option::isActive_){
		Cell *new_cell = new Cell(value, Row_access_option::columnIndex_, rowIndex);
		column.push_back(new_cell);
		Row_access_option::insert_cell(rowIndex, new_cell);
	} else {
		Cell *new_cell = new Cell(value, rowIndex);
		column.push_back(new_cell);
	}
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline void Vector_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::_update_cell(const Field_element_type &value, index rowIndex, index position)
{
	if constexpr (Row_access_option::isActive_){
		Cell *new_cell = new Cell(value, Row_access_option::columnIndex_, rowIndex);
		column_[position] = new_cell;
		Row_access_option::insert_cell(rowIndex, new_cell);
	} else {
		Cell *new_cell = new Cell(value, rowIndex);
		column_[position] = new_cell;
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // VECTOR_COLUMN_H
