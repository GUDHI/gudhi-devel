/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_HEAP_COLUMN_H
#define Z2_HEAP_COLUMN_H

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
	using Column_type = std::vector<Cell>;
	using iterator = typename Column_type::iterator;
	using const_iterator = typename Column_type::const_iterator;
	using reverse_iterator = typename Column_type::reverse_iterator;
	using const_reverse_iterator = typename Column_type::const_reverse_iterator;

	Z2_heap_column();
	template<class Container_type>
	Z2_heap_column(const Container_type& nonZeroRowIndices);
	template<class Container_type>
	Z2_heap_column(const Container_type& nonZeroRowIndices, dimension_type dimension);
	Z2_heap_column(const Z2_heap_column& column);
	Z2_heap_column(Z2_heap_column&& column) noexcept;

	std::vector<bool> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;
	reverse_iterator rbegin() noexcept;
	const_reverse_iterator rbegin() const noexcept;
	reverse_iterator rend() noexcept;
	const_reverse_iterator rend() const noexcept;

	Z2_heap_column& operator+=(Z2_heap_column const &column);
	friend Z2_heap_column operator+(Z2_heap_column column1, Z2_heap_column const& column2){
		column1 += column2;
		return column1;
	}

	Z2_heap_column& operator*=(unsigned int v);
	friend Z2_heap_column operator*(Z2_heap_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Z2_heap_column operator*(unsigned int const& v, Z2_heap_column column){
		column *= v;
		return column;
	}

	Z2_heap_column& operator=(Z2_heap_column other);

	friend void swap(Z2_heap_column& col1, Z2_heap_column& col2){
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
		std::swap(col1.insertsSinceLastPrune_, col2.insertsSinceLastPrune_);
	}

protected:
	int dim_;
	Column_type column_;
	unsigned int insertsSinceLastPrune_;

	void _prune();
	int _pop_pivot();
};

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column() : dim_(0), insertsSinceLastPrune_(0)
{}

template<class Column_pairing_option>
template<class Container_type>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column(const Container_type& nonZeroRowIndices)
	: dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
	  column_(nonZeroRowIndices.begin(), nonZeroRowIndices.end()),
	  insertsSinceLastPrune_(0)
{
	std::make_heap(column_.begin(), column_.end());
}

template<class Column_pairing_option>
template<class Container_type>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column(const Container_type& nonZeroRowIndices, dimension_type dimension)
	: dim_(dimension),
	  column_(nonZeroRowIndices.begin(), nonZeroRowIndices.end()),
	  insertsSinceLastPrune_(0)
{
	std::make_heap(column_.begin(), column_.end());
}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column(const Z2_heap_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_),
	  insertsSinceLastPrune_(column.insertsSinceLastPrune_)
{}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option>::Z2_heap_column(Z2_heap_column&& column) noexcept
	: Column_pairing_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_)),
	  insertsSinceLastPrune_(std::exchange(column.insertsSinceLastPrune_, 0))
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
inline typename Z2_heap_column<Column_pairing_option>::iterator
Z2_heap_column<Column_pairing_option>::begin() noexcept
{
	return column_.begin();
}

template<class Column_pairing_option>
inline typename Z2_heap_column<Column_pairing_option>::const_iterator
Z2_heap_column<Column_pairing_option>::begin() const noexcept
{
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
inline typename Z2_heap_column<Column_pairing_option>::reverse_iterator
Z2_heap_column<Column_pairing_option>::rbegin() noexcept
{
	return column_.rbegin();
}

template<class Column_pairing_option>
inline typename Z2_heap_column<Column_pairing_option>::const_reverse_iterator
Z2_heap_column<Column_pairing_option>::rbegin() const noexcept
{
	return column_.rbegin();
}

template<class Column_pairing_option>
inline typename Z2_heap_column<Column_pairing_option>::reverse_iterator
Z2_heap_column<Column_pairing_option>::rend() noexcept
{
	return column_.rend();
}

template<class Column_pairing_option>
inline typename Z2_heap_column<Column_pairing_option>::const_reverse_iterator
Z2_heap_column<Column_pairing_option>::rend() const noexcept
{
	return column_.rend();
}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option> &Z2_heap_column<Column_pairing_option>::operator+=(Z2_heap_column const &column)
{
	const Column_type& colToAdd = column.column_;
	const unsigned int size = colToAdd.size();

	if (size == 0) return *this;

	for (const Cell& v : colToAdd) {
		column_.push_back(v);
		std::push_heap(column_.begin(), column_.end());
	}
	insertsSinceLastPrune_ += size;

	if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

	return *this;
}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option> &Z2_heap_column<Column_pairing_option>::operator*=(unsigned int v)
{
	if (v % 2 == 0){
		column_.clear();
		insertsSinceLastPrune_ = 0;
	}

	return *this;
}

template<class Column_pairing_option>
inline Z2_heap_column<Column_pairing_option>& Z2_heap_column<Column_pairing_option>::operator=(Z2_heap_column other)
{
	std::swap(dim_, other.dim_);
	column_.swap(other.column_);
	std::swap(insertsSinceLastPrune_, other.insertsSinceLastPrune_);
	return *this;
}

template<class Column_pairing_option>
inline void Z2_heap_column<Column_pairing_option>::_prune()
{
	if (insertsSinceLastPrune_ == 0) return;

	Column_type tempCol;
	int pivot = _pop_pivot();
	while (pivot != -1) {
		tempCol.push_back(pivot);
		pivot = _pop_pivot();
	}
	column_.swap(tempCol);
	std::make_heap(column_.begin(), column_.end());

	insertsSinceLastPrune_ = 0;
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

	return pivot;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_HEAP_COLUMN_H
