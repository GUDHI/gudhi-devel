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

namespace Gudhi {
namespace persistence_matrix {

class Z2_set_column
{
public:
	Z2_set_column();
	Z2_set_column(boundary_type& boundary);
	Z2_set_column(Z2_set_column& column);
	Z2_set_column(Z2_set_column&& column) noexcept;

	void get_content(boundary_type& container);
	bool contains(unsigned int value) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	void clear();
	void clear(unsigned int value);
	void reorder(std::vector<index>& valueMap);
	void add(Z2_set_column& column);

	Z2_set_column& operator=(Z2_set_column other);

	friend void swap(Z2_set_column& col1, Z2_set_column& col2);

private:
	int dim_;
	std::set<unsigned int> column_;
};

inline Z2_set_column::Z2_set_column() : dim_(0)
{}

inline Z2_set_column::Z2_set_column(boundary_type &boundary)
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

inline void Z2_set_column::get_content(boundary_type &container)
{
	std::copy(column_.begin(), column_.end(), std::back_inserter(container));
}

inline bool Z2_set_column::contains(unsigned int value) const
{
	return column_.find(value) != column_.end();
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
	return *(column_.rbegin());
}

inline void Z2_set_column::clear()
{
	column_.clear();
}

inline void Z2_set_column::clear(unsigned int value)
{
	column_.erase(value);
}

inline void Z2_set_column::reorder(std::vector<index> &valueMap)
{
	std::set<unsigned int> newSet;
	for (const unsigned int& v : column_) newSet.insert(valueMap.at(v));
	column_.swap(newSet);
}

inline void Z2_set_column::add(Z2_set_column &column)
{
	for (const unsigned int& v : column.column_){
		if (column_.find(v) != column_.end())
			column_.erase(v);
		else
			column_.insert(v);
	}
}

inline Z2_set_column &Z2_set_column::operator=(Z2_set_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	return *this;
}

inline void swap(Z2_set_column& col1, Z2_set_column& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_SETCOLUMN_H
