/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_UNORDEREDSETCOLUMN_H
#define Z2_UNORDEREDSETCOLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>
#include <algorithm>

#include "../utilities/utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Column_pairing_option>
class Z2_unordered_set_column : public Column_pairing_option
{
public:
	using Cell = Z2_base_cell;
	using iterator = typename std::unordered_set<Cell>::iterator;
	using const_iterator = typename std::unordered_set<Cell>::const_iterator;

	Z2_unordered_set_column();
	template<class Boundary_type>
	Z2_unordered_set_column(Boundary_type& boundary);
	template<class Boundary_type>
	Z2_unordered_set_column(Boundary_type& boundary, dimension_type dimension);
	Z2_unordered_set_column(Z2_unordered_set_column& column);
	Z2_unordered_set_column(const Z2_unordered_set_column& column);
	Z2_unordered_set_column(Z2_unordered_set_column&& column) noexcept;

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

	Z2_unordered_set_column& operator+=(Z2_unordered_set_column const &column);
	friend Z2_unordered_set_column operator+(Z2_unordered_set_column column1,
											 Z2_unordered_set_column const& column2){
		column1 += column2;
		return column1;
	}

	Z2_unordered_set_column& operator=(Z2_unordered_set_column other);

	friend void swap(Z2_unordered_set_column& col1,
					 Z2_unordered_set_column& col2){
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
		std::swap(col1.pivotChanged_, col2.pivotChanged_);
		std::swap(col1.pivot_, col2.pivot_);
	}

private:
	int dim_;
	std::unordered_set<Cell> column_;
	bool pivotChanged_;
	int pivot_;
};

template<class Column_pairing_option>
inline Z2_unordered_set_column<Column_pairing_option>::Z2_unordered_set_column()
	: dim_(0), pivotChanged_(false), pivot_(-1)
{}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_unordered_set_column<Column_pairing_option>::Z2_unordered_set_column(Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.begin(), boundary.end()),
	  pivotChanged_(false),
	  pivot_(boundary.size() == 0 ? -1 : *(boundary.rbegin()))
{}

template<class Column_pairing_option>
template<class Boundary_type>
inline Z2_unordered_set_column<Column_pairing_option>::Z2_unordered_set_column(Boundary_type &boundary, dimension_type dimension)
	: dim_(dimension),
	  column_(boundary.begin(), boundary.end()),
	  pivotChanged_(false),
	  pivot_(boundary.size() == 0 ? -1 : *(boundary.rbegin()))
{}

template<class Column_pairing_option>
inline Z2_unordered_set_column<Column_pairing_option>::Z2_unordered_set_column(Z2_unordered_set_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_),
	  pivotChanged_(column.pivotChanged_),
	  pivot_(column.pivot_)
{}

template<class Column_pairing_option>
inline Z2_unordered_set_column<Column_pairing_option>::Z2_unordered_set_column(const Z2_unordered_set_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_),
	  column_(column.column_),
	  pivotChanged_(column.pivotChanged_),
	  pivot_(column.pivot_)
{}

template<class Column_pairing_option>
inline Z2_unordered_set_column<Column_pairing_option>::Z2_unordered_set_column(Z2_unordered_set_column &&column) noexcept
	: Column_pairing_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_)),
	  pivotChanged_(std::exchange(column.pivotChanged_, 0)),
	  pivot_(std::exchange(column.pivot_, 0))
{}

template<class Column_pairing_option>
inline std::vector<bool> Z2_unordered_set_column<Column_pairing_option>::get_content(unsigned int columnLength)
{
	std::vector<bool> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = 1;
	}
	return container;
}

template<class Column_pairing_option>
inline bool Z2_unordered_set_column<Column_pairing_option>::is_non_zero(index rowIndex) const
{
	return column_.find(rowIndex) != column_.end();
}

template<class Column_pairing_option>
inline bool Z2_unordered_set_column<Column_pairing_option>::is_empty()
{
	return column_.empty();
}

template<class Column_pairing_option>
inline dimension_type Z2_unordered_set_column<Column_pairing_option>::get_dimension() const
{
	return dim_;
}

template<class Column_pairing_option>
inline int Z2_unordered_set_column<Column_pairing_option>::get_pivot()
{
	if (pivotChanged_ && column_.size() == 0){
		pivot_ = -1;
		pivotChanged_ = false;
	} else if (pivotChanged_) {
		pivot_ = 0;
		for (const Cell& c : column_){
			if (static_cast<int>(c.get_row_index()) > pivot_)
				pivot_ = c.get_row_index();
		}
		pivotChanged_ = false;
	}

	return pivot_;
}

template<class Column_pairing_option>
inline void Z2_unordered_set_column<Column_pairing_option>::clear()
{
	column_.clear();
	pivot_ = -1;
	pivotChanged_ = false;
}

template<class Column_pairing_option>
inline void Z2_unordered_set_column<Column_pairing_option>::clear(index rowIndex)
{
	column_.erase(rowIndex);
	if (static_cast<int>(rowIndex) == pivot_) pivotChanged_ = true;
}

template<class Column_pairing_option>
inline void Z2_unordered_set_column<Column_pairing_option>::reorder(std::vector<index> &valueMap)
{
	std::unordered_set<Cell> newSet;
	for (const Cell& v : column_) newSet.insert(valueMap.at(v.get_row_index()));
	column_.swap(newSet);
	pivotChanged_ = true;
}

template<class Column_pairing_option>
inline typename Z2_unordered_set_column<Column_pairing_option>::iterator
Z2_unordered_set_column<Column_pairing_option>::begin() noexcept
{
	return column_.begin();
}

template<class Column_pairing_option>
inline typename Z2_unordered_set_column<Column_pairing_option>::const_iterator
Z2_unordered_set_column<Column_pairing_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Column_pairing_option>
inline typename Z2_unordered_set_column<Column_pairing_option>::iterator
Z2_unordered_set_column<Column_pairing_option>::end() noexcept
{
	return column_.end();
}

template<class Column_pairing_option>
inline typename Z2_unordered_set_column<Column_pairing_option>::const_iterator
Z2_unordered_set_column<Column_pairing_option>::end() const noexcept
{
	return column_.end();
}

template<class Column_pairing_option>
inline Z2_unordered_set_column<Column_pairing_option> &Z2_unordered_set_column<Column_pairing_option>::operator+=(Z2_unordered_set_column const &column)
{
	for (const Cell& v : column.column_){
		if (column_.find(v) != column_.end()){
			column_.erase(v);
			if (static_cast<int>(v.get_row_index()) == pivot_) pivotChanged_ = true;
		} else {
			column_.insert(v);
			if (static_cast<int>(v.get_row_index()) > pivot_){
				pivot_ = v.get_row_index();
				pivotChanged_ = false;
			}
		}
	}

	return *this;
}

template<class Column_pairing_option>
inline Z2_unordered_set_column<Column_pairing_option> &Z2_unordered_set_column<Column_pairing_option>::operator=(Z2_unordered_set_column other)
{
	std::swap(dim_, other.dim_);
	column_.swap(other.column_);
	std::swap(pivotChanged_, other.pivotChanged_);
	std::swap(pivot_, other.pivot_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_UNORDEREDSETCOLUMN_H
