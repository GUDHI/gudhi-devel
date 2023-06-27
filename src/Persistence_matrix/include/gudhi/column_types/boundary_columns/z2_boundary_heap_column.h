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

class Z2_heap_boundary_column : public Z2_heap_column
{
private:
	using Base = Z2_heap_column;
	using Base::operator+=;		//kinda ugly, so TODO: organize better

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

	std::vector<bool> get_content(int columnLength = -1) const;
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
		swap(static_cast<Z2_heap_column&>(col1),
			 static_cast<Z2_heap_column&>(col2));
		col1.erasedValues_.swap(col2.erasedValues_);
	}

private:
	std::unordered_set<unsigned int> erasedValues_;

	void _prune();
	int _pop_pivot();
};

inline Z2_heap_boundary_column::Z2_heap_boundary_column() : Base()
{}

template<class Boundary_type>
inline Z2_heap_boundary_column::Z2_heap_boundary_column(const Boundary_type& boundary)
	: Base(boundary)
{}

template<class Boundary_type>
inline Z2_heap_boundary_column::Z2_heap_boundary_column(const Boundary_type& boundary, dimension_type dimension)
	: Base(boundary, dimension)
{}

inline Z2_heap_boundary_column::Z2_heap_boundary_column(const Z2_heap_boundary_column &column)
	: Base(static_cast<const Base&>(column)),
	  erasedValues_(column.erasedValues_)
{}

inline Z2_heap_boundary_column::Z2_heap_boundary_column(Z2_heap_boundary_column&& column) noexcept
	: Base(std::move(static_cast<Base&&>(column))),
	  erasedValues_(std::move(column.erasedValues_))
{}

inline std::vector<bool> Z2_heap_boundary_column::get_content(int columnLength) const
{
	if (columnLength < 0) columnLength = Base::column_.front().get_row_index();

	std::vector<bool> container(columnLength, 0);
	for (auto it = Base::column_.begin(); it != Base::column_.end(); ++it){
		if (it->get_row_index() < static_cast<unsigned int>(columnLength) &&
				erasedValues_.find(it->get_row_index()) == erasedValues_.end())
			container[it->get_row_index()] = !container[it->get_row_index()];
	}
	return container;
}

inline bool Z2_heap_boundary_column::is_non_zero(index rowIndex) const
{
	if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

	return Base::is_non_zero(rowIndex);
}

inline bool Z2_heap_boundary_column::is_empty()
{
	int pivot = _pop_pivot();
	if (pivot != -1){
		Base::column_.push_back(pivot);
		std::push_heap(Base::column_.begin(), Base::column_.end());
		return false;
	}
	return true;
}

inline int Z2_heap_boundary_column::get_pivot()
{
	int pivot = _pop_pivot();
	if (pivot != -1){
		Base::column_.push_back(pivot);
		std::push_heap(Base::column_.begin(), Base::column_.end());
	}
	return pivot;
}

inline void Z2_heap_boundary_column::clear()
{
	Base::clear();
	erasedValues_.clear();
}

inline void Z2_heap_boundary_column::clear(index rowIndex)
{
	erasedValues_.insert(rowIndex);
}

template<class Map_type>
inline void Z2_heap_boundary_column::reorder(Map_type &valueMap)
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

inline typename Z2_heap_boundary_column::iterator
Z2_heap_boundary_column::begin() noexcept
{
	_prune();
	return Base::column_.begin();
}

inline Z2_heap_boundary_column &Z2_heap_boundary_column::operator+=(Z2_heap_boundary_column const &column)
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

inline Z2_heap_boundary_column &Z2_heap_boundary_column::operator*=(unsigned int v)
{
	if (v % 2 == 0){
		clear();
	}

	return *this;
}

inline Z2_heap_boundary_column& Z2_heap_boundary_column::operator=(Z2_heap_boundary_column other)
{
	Base::operator=(static_cast<Base&>(other));
	erasedValues_.swap(other.erasedValues_);
	return *this;
}

inline void Z2_heap_boundary_column::_prune()
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

inline int Z2_heap_boundary_column::_pop_pivot()
{
	unsigned int pivot = Base::_pop_pivot();

	if (erasedValues_.find(pivot) != erasedValues_.end())
		pivot = _pop_pivot();

	return pivot;
}

} //namespace persistence_matrix
} //namespace Gudhi

template<>
struct std::hash<Gudhi::persistence_matrix::Z2_heap_boundary_column>
{
	size_t operator()(const Gudhi::persistence_matrix::Z2_heap_boundary_column& column) const
	{
		std::size_t seed = 0;
		unsigned int i = 0;
		for (bool val : column.get_content()){
			seed ^= std::hash<unsigned int>()(i++ * val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

#endif // B_Z2_HEAP_COLUMN_H
