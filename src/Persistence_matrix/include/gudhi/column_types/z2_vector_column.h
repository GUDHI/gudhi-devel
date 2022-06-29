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

#include "../utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

class Z2_vector_column
{
public:
	using Cell = Z2_base_cell;

	Z2_vector_column();
	Z2_vector_column(boundary_type& boundary);
	Z2_vector_column(Z2_vector_column& column);
	Z2_vector_column(Z2_vector_column&& column) noexcept;

//	void get_content(boundary_type& container);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);
	void add(Z2_vector_column& column);

	Z2_vector_column& operator=(Z2_vector_column other);

	friend void swap(Z2_vector_column& col1, Z2_vector_column& col2);

private:
	int dim_;
	std::vector<Cell> column_;
	std::unordered_set<unsigned int> erasedValues_;

	void _cleanValues();
};

inline Z2_vector_column::Z2_vector_column() : dim_(0)
{}

inline Z2_vector_column::Z2_vector_column(boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.begin(), boundary.end())
{}

inline Z2_vector_column::Z2_vector_column(Z2_vector_column &column)
	: dim_(column.dim_),
	  column_(column.column_),
	  erasedValues_(column.erasedValues_)
{}

inline Z2_vector_column::Z2_vector_column(Z2_vector_column &&column) noexcept
	: dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_)),
	  erasedValues_(std::move(column.erasedValues_))
{}

//inline void Z2_vector_column::get_content(boundary_type &container)
//{
//	_cleanValues();
//	std::copy(column_.begin(), column_.end(), std::back_inserter(container));
//}

inline bool Z2_vector_column::is_non_zero(index rowIndex) const
{
	if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

	for (const Cell& v : column_){
		if (v.get_row_index() == rowIndex) return true;
	}
	return false;
}

inline bool Z2_vector_column::is_empty()
{
	_cleanValues();
	return column_.empty();
}

inline dimension_type Z2_vector_column::get_dimension() const
{
	return dim_;
}

inline int Z2_vector_column::get_pivot()
{
	auto it = erasedValues_.find(column_.back().get_row_index());
	while (!column_.empty() && it != erasedValues_.end()) {
		erasedValues_.erase(it);
		column_.pop_back();
		it = erasedValues_.find(column_.back().get_row_index());
	}

	if (column_.empty()) return -1;

	return column_.back().get_row_index();
}

inline void Z2_vector_column::clear()
{
	column_.clear();
	erasedValues_.clear();
}

inline void Z2_vector_column::clear(index rowIndex)
{
	erasedValues_.insert(rowIndex);
}

inline void Z2_vector_column::reorder(std::vector<index> &valueMap)
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

inline void Z2_vector_column::add(Z2_vector_column &column)
{
	if (column.is_empty()) return;
	if (column_.empty()){
		column._cleanValues();
		std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
		erasedValues_.clear();
		return;
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
}

inline Z2_vector_column &Z2_vector_column::operator=(Z2_vector_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	std::swap(erasedValues_, other.erasedValues_);
	return *this;
}

inline void Z2_vector_column::_cleanValues()
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

inline void swap(Z2_vector_column& col1, Z2_vector_column& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
	std::swap(col1.erasedValues_, col2.erasedValues_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_VECTORCOLUMN_H
