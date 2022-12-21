/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef B_Z2_HEAP_COLUMN_H
#define B_Z2_HEAP_COLUMN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "../../utilities/utilities.h"
#include "../z2_heap_column.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Column_pairing_option>
class Z2_heap_boundary_column : public Z2_heap_column<Column_pairing_option>
{
private:
	using Base = Z2_heap_column<Column_pairing_option>;

public:
	using Cell = typename Base::Cell;
	using Column_type = typename Base::Column_type;
	using iterator = typename Base::iterator;
	using const_iterator = typename Base::const_iterator;

	Z2_heap_boundary_column();
	template<class Boundary_type>
	Z2_heap_boundary_column(const Boundary_type& boundary);
	template<class Boundary_type>
	Z2_heap_boundary_column(const Boundary_type& boundary, dimension_type dimension);
	Z2_heap_boundary_column(const Z2_heap_boundary_column& column);
	Z2_heap_boundary_column(Z2_heap_boundary_column&& column) noexcept;

	std::vector<bool> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	template<class Map_type>
	void reorder(Map_type& valueMap);

	iterator begin() noexcept;

	Z2_heap_boundary_column& operator+=(Z2_heap_boundary_column const &column);
	friend Z2_heap_boundary_column operator+(Z2_heap_boundary_column column1, Z2_heap_boundary_column const& column2){
		column1 += column2;
		return column1;
	}
	Z2_heap_boundary_column& operator*=(unsigned int v);
	friend Z2_heap_boundary_column operator*(Z2_heap_boundary_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Z2_heap_boundary_column operator*(unsigned int const& v, Z2_heap_boundary_column column){
		column *= v;
		return column;
	}

	Z2_heap_boundary_column& operator=(Z2_heap_boundary_column other);

	friend void swap(Z2_heap_boundary_column& col1, Z2_heap_boundary_column& col2){
		swap(static_cast<Z2_heap_column<Column_pairing_option>&>(col1),
			 static_cast<Z2_heap_column<Column_pairing_option>&>(col2));
		col1.erasedValues_.swap(col2.erasedValues_);
	}

private:
	std::unordered_set<unsigned int> erasedValues_;

	void _prune();
	int _pop_pivot();
};

template<class Column_pairing_option>
inline Z2_heap_boundary_column<Column_pairing_option>::Z2_heap_boundary_column() : Base()
{}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_heap_boundary_column<Column_pairing_option>::Z2_heap_boundary_column(const Boundary_type& boundary)
	: Base(boundary)
{}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_heap_boundary_column<Column_pairing_option>::Z2_heap_boundary_column(const Boundary_type& boundary, dimension_type dimension)
	: Base(boundary, dimension)
{}

template<class Column_pairing_option>
inline Z2_heap_boundary_column<Column_pairing_option>::Z2_heap_boundary_column(const Z2_heap_boundary_column &column)
	: Base(static_cast<const Base&>(column)),
	  erasedValues_(column.erasedValues_)
{}

template<class Column_pairing_option>
inline Z2_heap_boundary_column<Column_pairing_option>::Z2_heap_boundary_column(Z2_heap_boundary_column&& column) noexcept
	: Base(std::move(static_cast<Base&&>(column))),
	  erasedValues_(std::move(column.erasedValues_))
{}

template<class Column_pairing_option>
inline std::vector<bool> Z2_heap_boundary_column<Column_pairing_option>::get_content(unsigned int columnLength)
{
	_prune();
	std::vector<bool> container(columnLength, 0);
	for (auto it = Base::column_.begin(); it != Base::column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = 1;
	}
	return container;
}

template<class Column_pairing_option>
inline bool Z2_heap_boundary_column<Column_pairing_option>::is_non_zero(index rowIndex) const
{
	if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

	return Base::is_non_zero(rowIndex);
}

template<class Column_pairing_option>
inline bool Z2_heap_boundary_column<Column_pairing_option>::is_empty()
{
	int pivot = _pop_pivot();
	if (pivot != -1){
		Base::column_.push_back(pivot);
		std::push_heap(Base::column_.begin(), Base::column_.end());
		return false;
	}
	return true;
}

template<class Column_pairing_option>
inline int Z2_heap_boundary_column<Column_pairing_option>::get_pivot()
{
	int pivot = _pop_pivot();
	if (pivot != -1){
		Base::column_.push_back(pivot);
		std::push_heap(Base::column_.begin(), Base::column_.end());
	}
	return pivot;
}

template<class Column_pairing_option>
inline void Z2_heap_boundary_column<Column_pairing_option>::clear()
{
	Base::column_.clear();
	Base::insertsSinceLastPrune_ = 0;
	erasedValues_.clear();
}

template<class Column_pairing_option>
inline void Z2_heap_boundary_column<Column_pairing_option>::clear(index rowIndex)
{
	erasedValues_.insert(rowIndex);
}

template<class Column_pairing_option>
template<class Map_type>
inline void Z2_heap_boundary_column<Column_pairing_option>::reorder(Map_type &valueMap)
{
	Column_type tempCol;
	int pivot = _pop_pivot();
	while (pivot != -1) {
		tempCol.push_back(valueMap[pivot]);
		pivot = _pop_pivot();
	}
	Base::column_.swap(tempCol);
	std::make_heap(Base::column_.begin(), Base::column_.end());

	Base::insertsSinceLastPrune_ = 0;
	erasedValues_.clear();
}

template<class Column_pairing_option>
inline typename Z2_heap_boundary_column<Column_pairing_option>::iterator
Z2_heap_boundary_column<Column_pairing_option>::begin() noexcept
{
	_prune();
	return Base::column_.begin();
}

template<class Column_pairing_option>
inline Z2_heap_boundary_column<Column_pairing_option> &Z2_heap_boundary_column<Column_pairing_option>::operator+=(Z2_heap_boundary_column const &column)
{
	const Column_type& colToAdd = column.column_;
	const unsigned int size = colToAdd.size();

	if (size == 0) return *this;

	for (const Cell& v : colToAdd) {
		if (column.erasedValues_.find(v.get_row_index()) == column.erasedValues_.end()){
			Base::column_.push_back(v);
			std::push_heap(Base::column_.begin(), Base::column_.end());
			erasedValues_.erase(v.get_row_index());
		}
	}
	Base::insertsSinceLastPrune_ += size;

	if (2 * Base::insertsSinceLastPrune_ > Base::column_.size()) _prune();

	return *this;
}

template<class Column_pairing_option>
inline Z2_heap_boundary_column<Column_pairing_option> &Z2_heap_boundary_column<Column_pairing_option>::operator*=(unsigned int v)
{
	if (v % 2 == 0){
		clear();
	}

	return *this;
}

template<class Column_pairing_option>
inline Z2_heap_boundary_column<Column_pairing_option>& Z2_heap_boundary_column<Column_pairing_option>::operator=(Z2_heap_boundary_column other)
{
	Base::operator=(static_cast<Base&>(other));
	erasedValues_.swap(other.erasedValues_);
	return *this;
}

template<class Column_pairing_option>
inline void Z2_heap_boundary_column<Column_pairing_option>::_prune()
{
	if (Base::insertsSinceLastPrune_ == 0 && erasedValues_.empty()) return;

	Column_type tempCol;
	int pivot = _pop_pivot();
	while (pivot != -1) {
		tempCol.push_back(pivot);
		pivot = _pop_pivot();
	}
	Base::column_.swap(tempCol);
	std::make_heap(Base::column_.begin(), Base::column_.end());

	Base::insertsSinceLastPrune_ = 0;
	erasedValues_.clear();
}

template<class Column_pairing_option>
inline int Z2_heap_boundary_column<Column_pairing_option>::_pop_pivot()
{
	unsigned int pivot = Base::_pop_pivot();

	if (erasedValues_.find(pivot) != erasedValues_.end())
		pivot = _pop_pivot();

	return pivot;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // B_Z2_HEAP_COLUMN_H
