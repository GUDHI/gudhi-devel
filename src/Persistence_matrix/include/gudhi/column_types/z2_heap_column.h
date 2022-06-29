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

#include "../utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

class Z2_heap_column
{
public:
	using Cell = Z2_base_cell;

	Z2_heap_column();
	Z2_heap_column(boundary_type& boundary);
	Z2_heap_column(Z2_heap_column& column);
	Z2_heap_column(Z2_heap_column&& column) noexcept;

//	void get_content(boundary_type& container);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);
	void add(Z2_heap_column& column);

	Z2_heap_column& operator=(Z2_heap_column other);

	friend void swap(Z2_heap_column& col1, Z2_heap_column& col2);

private:
	int dim_;
	std::vector<Cell> column_;
	unsigned int insertsSinceLastPrune_;
	std::unordered_set<unsigned int> erasedValues_;

	void _prune();
	int _pop_pivot();
};

inline Z2_heap_column::Z2_heap_column() : dim_(0), insertsSinceLastPrune_(0)
{}

inline Z2_heap_column::Z2_heap_column(boundary_type& boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.begin(), boundary.end()),
	  insertsSinceLastPrune_(0)
{
//	for (index i : boundary){
//		column_[i] = Cell(i);
//	}
	std::make_heap(column_.begin(), column_.end());
}

inline Z2_heap_column::Z2_heap_column(Z2_heap_column& column)
	: dim_(column.dim_),
	  column_(column.column_),
	  insertsSinceLastPrune_(column.insertsSinceLastPrune_),
	  erasedValues_(column.erasedValues_)
{}

inline Z2_heap_column::Z2_heap_column(Z2_heap_column&& column) noexcept
	: dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_)),
	  insertsSinceLastPrune_(std::exchange(column.insertsSinceLastPrune_, 0)),
	  erasedValues_(std::move(column.erasedValues_))
{}

//inline void Z2_heap_column::get_content(boundary_type &container)
//{
//	_prune();
//	container = column_;
//	std::sort_heap(container.begin(), container.end());
//}

inline bool Z2_heap_column::is_non_zero(index rowIndex) const
{
	if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

	unsigned int c = 0;

	for (const Cell& v : column_){
		if (v.get_row_index() == rowIndex) c++;
	}

	return c % 2 != 0;
}

inline bool Z2_heap_column::is_empty()
{
	int pivot = _pop_pivot();
	if (pivot != -1){
		column_.push_back(pivot);
		std::push_heap(column_.begin(), column_.end());
		return false;
	}
	return true;
}

inline dimension_type Z2_heap_column::get_dimension() const
{
	return dim_;
}

inline int Z2_heap_column::get_pivot()
{
	int pivot = _pop_pivot();
	if (pivot != -1){
		column_.push_back(pivot);
		std::push_heap(column_.begin(), column_.end());
	}
	return pivot;
}

inline void Z2_heap_column::clear()
{
	column_.clear();
	insertsSinceLastPrune_ = 0;
	erasedValues_.clear();
}

inline void Z2_heap_column::clear(index rowIndex)
{
	erasedValues_.insert(rowIndex);
}

inline void Z2_heap_column::reorder(std::vector<index> &valueMap)
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

inline void Z2_heap_column::add(Z2_heap_column &column)
{
	std::vector<Cell>& colToAdd = column.column_;
	const unsigned int size = colToAdd.size();

	if (size == 0) return;

	for (Cell& v : colToAdd) {
		if (column.erasedValues_.find(v.get_row_index()) == column.erasedValues_.end()){
			column_.push_back(v);
			std::push_heap(column_.begin(), column_.end());
			erasedValues_.erase(v.get_row_index());
		}
	}
	insertsSinceLastPrune_ += size;

	if (2 * insertsSinceLastPrune_ > column_.size()) _prune();
}

inline Z2_heap_column& Z2_heap_column::operator=(Z2_heap_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	std::swap(insertsSinceLastPrune_, other.insertsSinceLastPrune_);
	std::swap(erasedValues_, other.erasedValues_);
	return *this;
}

inline void Z2_heap_column::_prune()
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

inline int Z2_heap_column::_pop_pivot()
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

inline void swap(Z2_heap_column& col1, Z2_heap_column& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
	std::swap(col1.insertsSinceLastPrune_, col2.insertsSinceLastPrune_);
	std::swap(col1.erasedValues_, col2.erasedValues_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_HEAPCOLUMN_H
