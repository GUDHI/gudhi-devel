/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef VECTORCOLUMN_H
#define VECTORCOLUMN_H

#include <iostream>
#include <vector>
#include <unordered_set>

#include "../utilities/utilities.h"
#include "../utilities/Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Column_pairing_option>
class Vector_column : public Column_pairing_option
{
public:
	using Cell = Base_cell<Field_element_type>;
	using iterator = typename std::vector<Cell>::iterator;
	using const_iterator = typename std::vector<Cell>::const_iterator;

	Vector_column();
	template<class Boundary_type>
	Vector_column(const Boundary_type& boundary);
	template<class Boundary_type>
	Vector_column(const Boundary_type& boundary, dimension_type dimension);
//	Vector_column(Vector_column& column);
	Vector_column(const Vector_column& column);
	Vector_column(Vector_column&& column) noexcept;

	std::vector<Field_element_type> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	Field_element_type get_pivot_value();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;

	Vector_column& operator+=(Vector_column &column);
	friend Vector_column operator+(Vector_column column1, Vector_column& column2){
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
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
		col1.erasedValues_.swap(col2.erasedValues_);
	}

private:
	int dim_;
	std::vector<Cell> column_;
	std::unordered_set<unsigned int> erasedValues_;

	void _cleanValues();
};

template<class Field_element_type, class Column_pairing_option>
inline Vector_column<Field_element_type,Column_pairing_option>::Vector_column() : dim_(0)
{}

template<class Field_element_type, class Column_pairing_option>
template<class Boundary_type>
inline Vector_column<Field_element_type,Column_pairing_option>::Vector_column(const Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.size())
{
	for (unsigned int i = 0; i < boundary.size(); i++){
		column_[i] = Cell(boundary[i].second, boundary[i].first);
	}
}

template<class Field_element_type, class Column_pairing_option>
template<class Boundary_type>
inline Vector_column<Field_element_type,Column_pairing_option>::Vector_column(const Boundary_type &boundary, dimension_type dimension)
	: dim_(dimension),
	  column_(boundary.size())
{
	for (unsigned int i = 0; i < boundary.size(); i++){
		column_[i] = Cell(boundary[i].second, boundary[i].first);
	}
}

//template<class Field_element_type, class Column_pairing_option>
//inline Vector_column<Field_element_type,Column_pairing_option>::Vector_column(Vector_column &column)
//	: Column_pairing_option(column),
//	  dim_(column.dim_),
//	  column_(column.column_),
//	  erasedValues_(column.erasedValues_)
//{}

template<class Field_element_type, class Column_pairing_option>
inline Vector_column<Field_element_type,Column_pairing_option>::Vector_column(const Vector_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_),
	  erasedValues_(column.erasedValues_)
{}

template<class Field_element_type, class Column_pairing_option>
inline Vector_column<Field_element_type,Column_pairing_option>::Vector_column(Vector_column &&column) noexcept
	: Column_pairing_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_)),
	  erasedValues_(std::move(column.erasedValues_))
{}

template<class Field_element_type, class Column_pairing_option>
inline std::vector<Field_element_type> Vector_column<Field_element_type,Column_pairing_option>::get_content(unsigned int columnLength)
{
	_cleanValues();
	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = it->get_element();
	}
	return container;
}

template<class Field_element_type, class Column_pairing_option>
inline bool Vector_column<Field_element_type,Column_pairing_option>::is_non_zero(index rowIndex) const
{
	if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

	for (Cell v : column_){
		if (v.get_row_index() == rowIndex) return true;
	}
	return false;
}

template<class Field_element_type, class Column_pairing_option>
inline bool Vector_column<Field_element_type,Column_pairing_option>::is_empty()
{
	_cleanValues();
	return column_.empty();
}

template<class Field_element_type, class Column_pairing_option>
inline dimension_type Vector_column<Field_element_type,Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type, class Column_pairing_option>
inline int Vector_column<Field_element_type,Column_pairing_option>::get_pivot()
{
	while (!column_.empty() &&
		   erasedValues_.find(column_.back().get_row_index()) != erasedValues_.end()) {
		erasedValues_.erase(column_.back().get_row_index());
		column_.pop_back();
	}

	if (column_.empty()) return -1;

	return column_.back().get_row_index();
}

template<class Field_element_type, class Column_pairing_option>
inline Field_element_type Vector_column<Field_element_type,Column_pairing_option>::get_pivot_value()
{
	while (!column_.empty() &&
		   erasedValues_.find(column_.back().get_row_index()) != erasedValues_.end()) {
		erasedValues_.erase(column_.back().get_row_index());
		column_.pop_back();
	}

	if (column_.empty()) return 0;

	return column_.back().get_element();
}

template<class Field_element_type, class Column_pairing_option>
inline void Vector_column<Field_element_type,Column_pairing_option>::clear()
{
	column_.clear();
	erasedValues_.clear();
}

template<class Field_element_type, class Column_pairing_option>
inline void Vector_column<Field_element_type,Column_pairing_option>::clear(index rowIndex)
{
	erasedValues_.insert(rowIndex);
}

template<class Field_element_type, class Column_pairing_option>
inline void Vector_column<Field_element_type,Column_pairing_option>::reorder(std::vector<index> &valueMap)
{
	std::vector<Cell> newColumn;
	for (Cell& v : column_) {
		if (erasedValues_.find(v.get_row_index()) == erasedValues_.end())
			newColumn.push_back(Cell(v.get_element(), valueMap[v.get_row_index()]));
	}
	std::sort(newColumn.begin(), newColumn.end());
	erasedValues_.clear();
	column_.swap(newColumn);
}

template<class Field_element_type, class Column_pairing_option>
inline typename Vector_column<Field_element_type,Column_pairing_option>::iterator
Vector_column<Field_element_type,Column_pairing_option>::begin() noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Column_pairing_option>
inline typename Vector_column<Field_element_type,Column_pairing_option>::const_iterator
Vector_column<Field_element_type,Column_pairing_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Column_pairing_option>
inline typename Vector_column<Field_element_type,Column_pairing_option>::iterator
Vector_column<Field_element_type,Column_pairing_option>::end() noexcept
{
	return column_.end();
}

template<class Field_element_type, class Column_pairing_option>
inline typename Vector_column<Field_element_type,Column_pairing_option>::const_iterator
Vector_column<Field_element_type,Column_pairing_option>::end() const noexcept
{
	return column_.end();
}

template<class Field_element_type, class Column_pairing_option>
inline Vector_column<Field_element_type,Column_pairing_option> &Vector_column<Field_element_type,Column_pairing_option>::operator+=(Vector_column &column)
{
	if (column.is_empty()) return *this;

	if (column_.empty()){
		column._cleanValues();
		std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
		erasedValues_.clear();
		return *this;
	}

	std::vector<Cell> newColumn;

	typename std::vector<Cell>::iterator itToAdd = column.column_.begin();
	typename std::vector<Cell>::iterator itTarget = column_.begin();
	unsigned int curRowToAdd = itToAdd->get_row_index();
	unsigned int curRowTarget = itTarget->get_row_index();

	while (itToAdd != column.column_.end() && itTarget != column_.end())
	{
		while (itToAdd != column.column_.end() &&
			   column.erasedValues_.find(curRowToAdd) != column.erasedValues_.end()) {
			itToAdd++;
			curRowToAdd = itToAdd->get_row_index();
		}

		while (itTarget != column_.end() &&
			   erasedValues_.find(curRowTarget) != erasedValues_.end()) {
			itTarget++;
			curRowTarget = itTarget->get_row_index();
		}

		if (itToAdd != column.column_.end() && itTarget != column_.end()){
			if (curRowToAdd == curRowTarget){
				Field_element_type sum = itTarget->get_element() + itToAdd->get_element();
				if (sum != 0) newColumn.push_back(Cell(sum, curRowToAdd));
				itTarget++;
				itToAdd++;
			} else if (curRowToAdd < curRowTarget){
				newColumn.push_back(Cell(itToAdd->get_element(), curRowToAdd));
				itToAdd++;
			} else {
				newColumn.push_back(Cell(itTarget->get_element(), curRowTarget));
				itTarget++;
			}
		}

		if (itToAdd != column.column_.end()) curRowToAdd = itToAdd->get_row_index();
		if (itTarget != column_.end()) curRowTarget = itTarget->get_row_index();
	}

	while (itToAdd != column.column_.end()){
		while (itToAdd != column.column_.end() &&
			   column.erasedValues_.find(curRowToAdd) != column.erasedValues_.end()) {
			itToAdd++;
			curRowToAdd = itToAdd->get_row_index();
		}

		if (itToAdd != column.column_.end()){
			newColumn.push_back(Cell(itToAdd->get_element(), curRowToAdd));
			itToAdd++;
		}

		if (itToAdd != column.column_.end()) curRowToAdd = itToAdd->get_row_index();
	}

	while (itTarget != column_.end()){
		while (itTarget != column_.end() &&
			   erasedValues_.find(curRowTarget) != erasedValues_.end()) {
			itTarget++;
			curRowTarget = itTarget->get_row_index();
		}

		if (itTarget != column_.end()){
			newColumn.push_back(Cell(itTarget->get_element(), curRowTarget));
			itTarget++;
		}

		if (itTarget != column_.end()) curRowTarget = itTarget->get_row_index();
	}

	column_.swap(newColumn);
	erasedValues_.clear();

	return *this;
}

template<class Field_element_type, class Column_pairing_option>
inline Vector_column<Field_element_type,Column_pairing_option> &Vector_column<Field_element_type,Column_pairing_option>::operator*=(unsigned int v)
{
	v %= Field_element_type::get_characteristic();

	if (v == 0) {
		column_.clear();
		erasedValues_.clear();
		return *this;
	}

	for (Cell& cell : column_){
		cell.get_element() *= v;
	}

	return *this;
}

template<class Field_element_type, class Column_pairing_option>
inline Vector_column<Field_element_type,Column_pairing_option> &Vector_column<Field_element_type,Column_pairing_option>::operator=(Vector_column other)
{
	std::swap(dim_, other.dim_);
	column_.swap(other.column_);
	erasedValues_.swap(other.erasedValues_);
	return *this;
}

template<class Field_element_type, class Column_pairing_option>
inline void Vector_column<Field_element_type,Column_pairing_option>::_cleanValues()
{
	if (erasedValues_.empty()) return;

	std::vector<Cell> newColumn;
	for (Cell& v : column_){
		if (erasedValues_.find(v.get_row_index()) == erasedValues_.end())
			newColumn.push_back(v);
	}
	erasedValues_.clear();
	column_.swap(newColumn);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // VECTORCOLUMN_H
