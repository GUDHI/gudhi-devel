/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_SETCOLUMN_H
#define Z2_SETCOLUMN_H

#include <iostream>
#include <list>
#include <set>

#include "../utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

class Z2_set_column
{
public:
	using Cell = Z2_base_cell;

	Z2_set_column();
	template<class Boundary_type>
	Z2_set_column(Boundary_type& boundary);
	Z2_set_column(Z2_set_column& column);
	Z2_set_column(Z2_set_column&& column) noexcept;

	std::vector<bool> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);

	Z2_set_column& operator+=(Z2_set_column const &column);
	friend Z2_set_column operator+(Z2_set_column column1, Z2_set_column const& column2);

	Z2_set_column& operator=(Z2_set_column other);

	friend void swap(Z2_set_column& col1, Z2_set_column& col2);

private:
	int dim_;
	std::set<Cell> column_;
};

inline Z2_set_column::Z2_set_column() : dim_(0)
{}

template<class Boundary_type>
inline Z2_set_column::Z2_set_column(Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.begin(), boundary.end())
{}

inline Z2_set_column::Z2_set_column(Z2_set_column &column)
	: dim_(column.dim_),
	  column_(column.column_)
{}

inline Z2_set_column::Z2_set_column(Z2_set_column &&column) noexcept
	: dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

inline std::vector<bool> Z2_set_column::get_content(unsigned int columnLength)
{
	std::vector<bool> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = 1;
	}
	return container;
}

inline bool Z2_set_column::is_non_zero(index rowIndex) const
{
	return column_.find(rowIndex) != column_.end();
}

inline bool Z2_set_column::is_empty()
{
	return column_.empty();
}

inline dimension_type Z2_set_column::get_dimension() const
{
	return dim_;
}

inline int Z2_set_column::get_pivot()
{
	if (column_.empty()) return -1;
	return column_.rbegin()->get_row_index();
}

inline void Z2_set_column::clear()
{
	column_.clear();
}

inline void Z2_set_column::clear(index rowIndex)
{
	column_.erase(rowIndex);
}

inline void Z2_set_column::reorder(std::vector<index> &valueMap)
{
	std::set<Cell> newSet;
	for (const Cell& v : column_) newSet.insert(valueMap.at(v.get_row_index()));
	column_.swap(newSet);
}

inline Z2_set_column &Z2_set_column::operator+=(Z2_set_column const &column)
{
	for (const Cell& v : column.column_){
		if (column_.find(v) != column_.end())
			column_.erase(v);
		else
			column_.insert(v);
	}

	return *this;
}

inline Z2_set_column &Z2_set_column::operator=(Z2_set_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	return *this;
}

Z2_set_column operator+(Z2_set_column column1, Z2_set_column const& column2)
{
	column1 += column2;
	return column1;
}

inline void swap(Z2_set_column& col1, Z2_set_column& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_SETCOLUMN_H
