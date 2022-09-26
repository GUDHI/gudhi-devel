/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_HEAPCOLUMN_H
#define Z2_HEAPCOLUMN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "../utilities/utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Column_pairing_option>
class Z2_heap_column : public Column_pairing_option
{
public:
	using Cell = Z2_base_cell;
	using iterator = typename std::vector<Cell>::iterator;
	using const_iterator = typename std::vector<Cell>::const_iterator;

	Z2_heap_column();
	template<class Boundary_type>
	Z2_heap_column(Boundary_type& boundary);
	template<class Boundary_type>
	Z2_heap_column(Boundary_type& boundary, dimension_type dimension);
	Z2_heap_column(Z2_heap_column& column);
	Z2_heap_column(const Z2_heap_column& column);
	Z2_heap_column(Z2_heap_column&& column) noexcept;

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

	Z2_heap_column& operator+=(Z2_heap_column const &column);
	template<class Friend_column_pairing_option>
	friend Z2_heap_column<Friend_column_pairing_option> operator+(
			Z2_heap_column<Friend_column_pairing_option> column1,
			Z2_heap_column<Friend_column_pairing_option> const& column2);

	Z2_heap_column& operator=(Z2_heap_column other);

	template<class Friend_column_pairing_option>
	friend void swap(
			Z2_heap_column<Friend_column_pairing_option>& col1,
			Z2_heap_column<Friend_column_pairing_option>& col2);

private:
	int dim_;
	std::vector<Cell> column_;
	unsigned int insertsSinceLastPrune_;
	std::unordered_set<unsigned int> erasedValues_;

	void _prune();
	int _pop_pivot();
};

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column() : dim_(0), insertsSinceLastPrune_(0)
{}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column(Boundary_type& boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.begin(), boundary.end()),
	  insertsSinceLastPrune_(0)
{
	std::make_heap(column_.begin(), column_.end());
}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column(Boundary_type& boundary, dimension_type dimension)
	: dim_(dimension),
	  column_(boundary.begin(), boundary.end()),
	  insertsSinceLastPrune_(0)
{
	std::make_heap(column_.begin(), column_.end());
}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column(Z2_heap_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_),
	  insertsSinceLastPrune_(column.insertsSinceLastPrune_),
	  erasedValues_(column.erasedValues_)
{}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column(const Z2_heap_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_),
	  insertsSinceLastPrune_(column.insertsSinceLastPrune_),
	  erasedValues_(column.erasedValues_)
{}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column(Z2_heap_column&& column) noexcept
	: Column_pairing_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_)),
	  insertsSinceLastPrune_(std::exchange(column.insertsSinceLastPrune_, 0)),
	  erasedValues_(std::move(column.erasedValues_))
{}

template<class Column_pairing_option>
inline std::vector<bool> Z2_heap_column<Column_pairing_option>::get_content(unsigned int columnLength)
{
	_prune();
	std::vector<bool> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = 1;
	}
	return container;
}

template<class Column_pairing_option>
inline bool Z2_heap_column<Column_pairing_option>::is_non_zero(index rowIndex) const
{
	if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

	unsigned int c = 0;

	for (const Cell& v : column_){
		if (v.get_row_index() == rowIndex) c++;
	}

	return c % 2 != 0;
}

template<class Column_pairing_option>
inline bool Z2_heap_column<Column_pairing_option>::is_empty()
{
	int pivot = _pop_pivot();
	if (pivot != -1){
		column_.push_back(pivot);
		std::push_heap(column_.begin(), column_.end());
		return false;
	}
	return true;
}

template<class Column_pairing_option>
inline dimension_type Z2_heap_column<Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Column_pairing_option>
inline int Z2_heap_column<Column_pairing_option>::get_pivot()
{
	int pivot = _pop_pivot();
	if (pivot != -1){
		column_.push_back(pivot);
		std::push_heap(column_.begin(), column_.end());
	}
	return pivot;
}

template<class Column_pairing_option>
inline void Z2_heap_column<Column_pairing_option>::clear()
{
	column_.clear();
	insertsSinceLastPrune_ = 0;
	erasedValues_.clear();
}

template<class Column_pairing_option>
inline void Z2_heap_column<Column_pairing_option>::clear(index rowIndex)
{
	erasedValues_.insert(rowIndex);
}

template<class Column_pairing_option>
inline void Z2_heap_column<Column_pairing_option>::reorder(std::vector<index> &valueMap)
{
	std::vector<Cell> tempCol;
	int pivot = _pop_pivot();
	while (pivot != -1) {
		tempCol.push_back(valueMap.at(pivot));
		pivot = _pop_pivot();
	}
	column_.swap(tempCol);
	std::make_heap(column_.begin(), column_.end());

	insertsSinceLastPrune_ = 0;
	erasedValues_.clear();
}

template<class Column_pairing_option>
inline typename Z2_heap_column<Column_pairing_option>::iterator
Z2_heap_column<Column_pairing_option>::begin() noexcept
{
	_prune();
	return column_.begin();
}

template<class Column_pairing_option>
inline typename Z2_heap_column<Column_pairing_option>::const_iterator
Z2_heap_column<Column_pairing_option>::begin() const noexcept
{
	_prune();
	return column_.begin();
}

template<class Column_pairing_option>
inline typename Z2_heap_column<Column_pairing_option>::iterator
Z2_heap_column<Column_pairing_option>::end() noexcept
{
	return column_.end();
}

template<class Column_pairing_option>
inline typename Z2_heap_column<Column_pairing_option>::const_iterator
Z2_heap_column<Column_pairing_option>::end() const noexcept
{
	return column_.end();
}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option> &Z2_heap_column<Column_pairing_option>::operator+=(Z2_heap_column const &column)
{
	const std::vector<Cell>& colToAdd = column.column_;
	const unsigned int size = colToAdd.size();

	if (size == 0) return *this;

	for (const Cell& v : colToAdd) {
		if (column.erasedValues_.find(v.get_row_index()) == column.erasedValues_.end()){
			column_.push_back(v);
			std::push_heap(column_.begin(), column_.end());
			erasedValues_.erase(v.get_row_index());
		}
	}
	insertsSinceLastPrune_ += size;

	if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

	return *this;
}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option>& Z2_heap_column<Column_pairing_option>::operator=(Z2_heap_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	std::swap(insertsSinceLastPrune_, other.insertsSinceLastPrune_);
	std::swap(erasedValues_, other.erasedValues_);
	return *this;
}

template<class Column_pairing_option>
inline void Z2_heap_column<Column_pairing_option>::_prune()
{
	if (insertsSinceLastPrune_ == 0 && erasedValues_.empty()) return;

	std::vector<Cell> tempCol;
	int pivot = _pop_pivot();
	while (pivot != -1) {
		tempCol.push_back(pivot);
		pivot = _pop_pivot();
	}
	column_.swap(tempCol);
	std::make_heap(column_.begin(), column_.end());

	insertsSinceLastPrune_ = 0;
	erasedValues_.clear();
}

template<class Column_pairing_option>
inline int Z2_heap_column<Column_pairing_option>::_pop_pivot()
{
	if (column_.empty()) {
		return -1;
	}

	unsigned int pivot = column_.front().get_row_index();
	std::pop_heap(column_.begin(), column_.end());
	column_.pop_back();
	while (!column_.empty() && column_.front().get_row_index() == pivot)
	{
		std::pop_heap(column_.begin(), column_.end());
		column_.pop_back();

		if (column_.empty()) {
			return -1;
		}
		pivot = column_.front().get_row_index();
		std::pop_heap(column_.begin(), column_.end());
		column_.pop_back();
	}

	if (erasedValues_.find(pivot) != erasedValues_.end())
		pivot = _pop_pivot();

	return pivot;
}

template <class Friend_column_pairing_option>
Z2_heap_column<Friend_column_pairing_option> operator+(
		Z2_heap_column<Friend_column_pairing_option> column1,
		Z2_heap_column<Friend_column_pairing_option> const& column2)
{
	column1 += column2;
	return column1;
}

template<class Friend_column_pairing_option>
inline void swap(Z2_heap_column<Friend_column_pairing_option>& col1,
				 Z2_heap_column<Friend_column_pairing_option>& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
	std::swap(col1.insertsSinceLastPrune_, col2.insertsSinceLastPrune_);
	std::swap(col1.erasedValues_, col2.erasedValues_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_HEAPCOLUMN_H
