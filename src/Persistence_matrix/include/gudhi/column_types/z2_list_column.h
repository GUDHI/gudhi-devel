/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_LISTCOLUMN_H
#define Z2_LISTCOLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>

#include "../utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

class Z2_list_column
{
public:
	using Cell = Z2_base_cell;

	Z2_list_column();
	template<class Boundary_type>
	Z2_list_column(Boundary_type& boundary);
	Z2_list_column(Z2_list_column& column);
	Z2_list_column(Z2_list_column&& column) noexcept;

	std::vector<bool> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	dimension_type get_dimension() const;
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);

	Z2_list_column& operator+=(Z2_list_column const &column);
	friend Z2_list_column operator+(Z2_list_column column1, Z2_list_column const& column2);

	Z2_list_column& operator=(Z2_list_column other);

	friend void swap(Z2_list_column& col1, Z2_list_column& col2);

private:
	int dim_;
	std::list<Cell> column_;
};

inline Z2_list_column::Z2_list_column() : dim_(0)
{}

template<class Boundary_type>
inline Z2_list_column::Z2_list_column(Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.begin(), boundary.end())
{}

inline Z2_list_column::Z2_list_column(Z2_list_column &column)
	: dim_(column.dim_),
	  column_(column.column_)
{}

inline Z2_list_column::Z2_list_column(Z2_list_column &&column) noexcept
	: dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

inline std::vector<bool> Z2_list_column::get_content(unsigned int columnLength)
{
	std::vector<bool> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = 1;
	}
	return container;
}

inline bool Z2_list_column::is_non_zero(index rowIndex) const
{
	for (const Cell& v : column_){
		if (v.get_row_index() == rowIndex) return true;
	}
	return false;
}

inline bool Z2_list_column::is_empty() const
{
	return column_.empty();
}

inline dimension_type Z2_list_column::get_dimension() const
{
	return dim_;
}

inline int Z2_list_column::get_pivot()
{
	if (column_.empty()) return -1;

	return column_.back().get_row_index();
}

inline void Z2_list_column::clear()
{
	column_.clear();
}

inline void Z2_list_column::clear(index rowIndex)
{
	auto it = column_.begin();
	while (it != column_.end() && it->get_row_index() != rowIndex) it++;
	if (it != column_.end()) column_.erase(it);
}

inline void Z2_list_column::reorder(std::vector<index> &valueMap)
{
	std::list<Cell>::iterator it = column_.begin();
	while (it != column_.end()) {
		it->setRowIndex(valueMap.at(it->get_row_index()));
		it++;
	}
	column_.sort();
}

inline Z2_list_column &Z2_list_column::operator+=(Z2_list_column const &column)
{
	if (column.is_empty()) return *this;
	if (column_.empty()){
		std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
		return *this;
	}

	std::list<Cell>::const_iterator itToAdd = column.column_.begin();
	std::list<Cell>::iterator itTarget = column_.begin();
	unsigned int valToAdd = itToAdd->get_row_index();
	unsigned int valTarget = itTarget->get_row_index();

	while (itToAdd != column.column_.end() && itTarget != column_.end())
	{
		if (itToAdd != column.column_.end() && itTarget != column_.end()){
			if (valToAdd == valTarget){
				column_.erase(itTarget++);
				itToAdd++;
			} else if (valToAdd < valTarget){
				column_.insert(itTarget, valToAdd);
				itToAdd++;
			} else {
				itTarget++;
			}
		}

		valToAdd = itToAdd->get_row_index();
		valTarget = itTarget->get_row_index();
	}

	while (itToAdd != column.column_.end()){
		valToAdd = itToAdd->get_row_index();
		if (itToAdd != column.column_.end()){
			column_.push_back(valToAdd);
			itToAdd++;
		}
	}

	return *this;
}

inline Z2_list_column &Z2_list_column::operator=(Z2_list_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	return *this;
}

Z2_list_column operator+(Z2_list_column column1, Z2_list_column const& column2)
{
	column1 += column2;
	return column1;
}

inline void swap(Z2_list_column& col1, Z2_list_column& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_LISTCOLUMN_H
