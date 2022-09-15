/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_VECTORCOLUMN_H
#define Z2_VECTORCOLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>

#include "../utilities/utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Column_pairing_option>
class Z2_vector_column : public Column_pairing_option
{
public:
	using Cell = Z2_base_cell;
	using iterator = typename std::vector<Cell>::iterator;
	using const_iterator = typename std::vector<Cell>::const_iterator;

	Z2_vector_column();
	template<class Boundary_type>
	Z2_vector_column(Boundary_type& boundary);
	template<class Boundary_type>
	Z2_vector_column(Boundary_type& boundary, dimension_type dimension);
	Z2_vector_column(Z2_vector_column& column);
	Z2_vector_column(const Z2_vector_column& column);
	Z2_vector_column(Z2_vector_column&& column) noexcept;

	std::vector<bool> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;

	Z2_vector_column& operator+=(Z2_vector_column &column);
	template<class Friend_column_pairing_option>
	friend Z2_vector_column<Friend_column_pairing_option> operator+(
			Z2_vector_column<Friend_column_pairing_option> column1,
			Z2_vector_column<Friend_column_pairing_option>& column2);

	Z2_vector_column& operator=(Z2_vector_column other);

	template<class Friend_column_pairing_option>
	friend void swap(Z2_vector_column<Friend_column_pairing_option>& col1,
					 Z2_vector_column<Friend_column_pairing_option>& col2);

private:
	int dim_;
	std::vector<Cell> column_;
	std::unordered_set<unsigned int> erasedValues_;

	void _cleanValues();
};

template<class Column_pairing_option>
inline Z2_vector_column<Column_pairing_option>::Z2_vector_column() : dim_(0)
{}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_vector_column<Column_pairing_option>::Z2_vector_column(Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.begin(), boundary.end())
{}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_vector_column<Column_pairing_option>::Z2_vector_column(Boundary_type &boundary, dimension_type dimension)
	: dim_(dimension),
	  column_(boundary.begin(), boundary.end())
{}

template<class Column_pairing_option>
inline Z2_vector_column<Column_pairing_option>::Z2_vector_column(Z2_vector_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_),
	  erasedValues_(column.erasedValues_)
{}

template<class Column_pairing_option>
inline Z2_vector_column<Column_pairing_option>::Z2_vector_column(const Z2_vector_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_),
	  erasedValues_(column.erasedValues_)
{}

template<class Column_pairing_option>
inline Z2_vector_column<Column_pairing_option>::Z2_vector_column(Z2_vector_column &&column) noexcept
	: Column_pairing_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_)),
	  erasedValues_(std::move(column.erasedValues_))
{}

template<class Column_pairing_option>
inline std::vector<bool> Z2_vector_column<Column_pairing_option>::get_content(unsigned int columnLength)
{
	_cleanValues();
	std::vector<bool> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = 1;
	}
	return container;
}

template<class Column_pairing_option>
inline bool Z2_vector_column<Column_pairing_option>::is_non_zero(index rowIndex) const
{
	if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

	for (const Cell& v : column_){
		if (v.get_row_index() == rowIndex) return true;
	}
	return false;
}

template<class Column_pairing_option>
inline bool Z2_vector_column<Column_pairing_option>::is_empty()
{
	_cleanValues();
	return column_.empty();
}

template<class Column_pairing_option>
inline dimension_type Z2_vector_column<Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Column_pairing_option>
inline int Z2_vector_column<Column_pairing_option>::get_pivot()
{
	if (column_.empty()) return -1;

	auto it = erasedValues_.find(column_.back().get_row_index());
	while (!column_.empty() && it != erasedValues_.end()) {
		erasedValues_.erase(it);
		column_.pop_back();
		it = erasedValues_.find(column_.back().get_row_index());
	}

	if (column_.empty()) return -1;

	return column_.back().get_row_index();
}

template<class Column_pairing_option>
inline void Z2_vector_column<Column_pairing_option>::clear()
{
	column_.clear();
	erasedValues_.clear();
}

template<class Column_pairing_option>
inline void Z2_vector_column<Column_pairing_option>::clear(index rowIndex)
{
	erasedValues_.insert(rowIndex);
}

template<class Column_pairing_option>
inline void Z2_vector_column<Column_pairing_option>::reorder(std::vector<index> &valueMap)
{
	std::vector<Cell> newColumn;
	for (const Cell& v : column_) {
		if (erasedValues_.find(v.get_row_index()) == erasedValues_.end())
			newColumn.push_back(valueMap.at(v.get_row_index()));
	}
	std::sort(newColumn.begin(), newColumn.end());
	erasedValues_.clear();
	column_.swap(newColumn);
}

template<class Column_pairing_option>
inline typename Z2_vector_column<Column_pairing_option>::iterator
Z2_vector_column<Column_pairing_option>::begin() noexcept
{
	return column_.begin();
}

template<class Column_pairing_option>
inline typename Z2_vector_column<Column_pairing_option>::const_iterator
Z2_vector_column<Column_pairing_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Column_pairing_option>
inline typename Z2_vector_column<Column_pairing_option>::iterator
Z2_vector_column<Column_pairing_option>::end() noexcept
{
	return column_.end();
}

template<class Column_pairing_option>
inline typename Z2_vector_column<Column_pairing_option>::const_iterator
Z2_vector_column<Column_pairing_option>::end() const noexcept
{
	return column_.end();
}

template<class Column_pairing_option>
inline Z2_vector_column<Column_pairing_option> &Z2_vector_column<Column_pairing_option>::operator+=(Z2_vector_column &column)
{
	if (column.is_empty()) return *this;
	if (column_.empty()){
		column._cleanValues();
		std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
		erasedValues_.clear();
		return *this;
	}

	std::vector<Cell> newColumn;

	std::vector<Cell>::iterator itToAdd = column.column_.begin();
	std::vector<Cell>::iterator itTarget = column_.begin();
	unsigned int valToAdd = itToAdd->get_row_index();
	unsigned int valTarget = itTarget->get_row_index();

	while (itToAdd != column.column_.end() && itTarget != column_.end())
	{
		while (itToAdd != column.column_.end() &&
			   column.erasedValues_.find(valToAdd) != column.erasedValues_.end()) {
			itToAdd++;
			valToAdd = itToAdd->get_row_index();
		}

		while (itTarget != column_.end() &&
			   erasedValues_.find(valTarget) != erasedValues_.end()) {
			itTarget++;
			valTarget = itTarget->get_row_index();
		}

		if (itToAdd != column.column_.end() && itTarget != column_.end()){
			if (valToAdd == valTarget){
				itTarget++;
				itToAdd++;
			} else if (valToAdd < valTarget){
				newColumn.push_back(valToAdd);
				itToAdd++;
			} else {
				newColumn.push_back(valTarget);
				itTarget++;
			}
		}

		valToAdd = itToAdd->get_row_index();
		valTarget = itTarget->get_row_index();
	}

	while (itToAdd != column.column_.end()){
		while (itToAdd != column.column_.end() &&
			   column.erasedValues_.find(valToAdd) != column.erasedValues_.end()) {
			itToAdd++;
			valToAdd = itToAdd->get_row_index();
		}

		if (itToAdd != column.column_.end()){
			newColumn.push_back(*itToAdd);
			itToAdd++;
		}
	}

	while (itTarget != column_.end()){
		while (itTarget != column_.end() &&
			   erasedValues_.find(valTarget) != erasedValues_.end()) {
			itTarget++;
			valTarget = itTarget->get_row_index();
		}

		if (itTarget != column_.end()){
			newColumn.push_back(*itTarget);
			itTarget++;
		}
	}

	column_.swap(newColumn);
	erasedValues_.clear();

	return *this;
}

template<class Column_pairing_option>
inline Z2_vector_column<Column_pairing_option> &Z2_vector_column<Column_pairing_option>::operator=(Z2_vector_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	std::swap(erasedValues_, other.erasedValues_);
	return *this;
}

template<class Column_pairing_option>
inline void Z2_vector_column<Column_pairing_option>::_cleanValues()
{
	if (erasedValues_.empty()) return;

	std::vector<Cell> newColumn;
	for (const Cell& v : column_){
		if (erasedValues_.find(v.get_row_index()) == erasedValues_.end())
			newColumn.push_back(v);
	}
	erasedValues_.clear();
	column_.swap(newColumn);
}

template<class Friend_column_pairing_option>
Z2_vector_column<Friend_column_pairing_option> operator+(
		Z2_vector_column<Friend_column_pairing_option> column1,
		Z2_vector_column<Friend_column_pairing_option>& column2)
{
	column1 += column2;
	return column1;
}

template<class Friend_column_pairing_option>
inline void swap(Z2_vector_column<Friend_column_pairing_option>& col1,
				 Z2_vector_column<Friend_column_pairing_option>& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
	std::swap(col1.erasedValues_, col2.erasedValues_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_VECTORCOLUMN_H
