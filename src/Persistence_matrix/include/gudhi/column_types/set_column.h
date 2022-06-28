/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SETCOLUMN_H
#define SETCOLUMN_H

#include <iostream>
#include <list>
#include <set>

#include "../utilities.h"
#include "../Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type = Zp_field_element<11> >
class Set_column
{
public:
	using Cell = Base_cell<Field_element_type>;

	Set_column();
	Set_column(std::vector<index>& rowIndices, std::vector<unsigned int>& values);
	Set_column(Set_column& column);
	Set_column(Set_column&& column) noexcept;

//	void get_content(boundary_type& container);
	bool isNonZero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);
	void add(Set_column& column);

	Set_column& operator=(Set_column other);

	template<class Friend_field_element_type>
	friend void swap(Set_column<Friend_field_element_type>& col1,
					 Set_column<Friend_field_element_type>& col2);

private:
	int dim_;
	std::set<Cell> column_;
};

template<class Field_element_type>
inline Set_column<Field_element_type>::Set_column() : dim_(0)
{}

template<class Field_element_type>
inline Set_column<Field_element_type>::Set_column(
		std::vector<index> &rowIndices, std::vector<unsigned int> &values)
	: dim_(rowIndices.size() == 0 ? 0 : rowIndices.size() - 1),
	  column_(rowIndices.size())
{
	for (unsigned int i = 0; i < rowIndices.size(); ++i){
		column_.insert(Cell(values.at(i), rowIndices[i]));
	}
}

template<class Field_element_type>
inline Set_column<Field_element_type>::Set_column(Set_column &column)
	: dim_(column.dim_),
	  column_(column.column_)
{}

template<class Field_element_type>
inline Set_column<Field_element_type>::Set_column(Set_column &&column) noexcept
	: dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

//template<class Field_element_type>
//inline void Set_column<Field_element_type>::get_content(boundary_type &container)
//{
//	std::copy(column_.begin(), column_.end(), std::back_inserter(container));
//}

template<class Field_element_type>
inline bool Set_column<Field_element_type>::isNonZero(index rowIndex) const
{
	return column_.find(Cell(0, rowIndex)) != column_.end();
}

template<class Field_element_type>
inline bool Set_column<Field_element_type>::is_empty()
{
	return column_.empty();
}

template<class Field_element_type>
inline dimension_type Set_column<Field_element_type>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type>
inline int Set_column<Field_element_type>::get_pivot()
{
	if (column_.empty()) return -1;
	return column_.rbegin()->get_row_index();
}

template<class Field_element_type>
inline void Set_column<Field_element_type>::clear()
{
	column_.clear();
}

template<class Field_element_type>
inline void Set_column<Field_element_type>::clear(index rowIndex)
{
	column_.erase(Cell(0, rowIndex));
}

template<class Field_element_type>
inline void Set_column<Field_element_type>::reorder(std::vector<index> &valueMap)
{
	std::set<Cell> newSet;
	for (const Cell& v : column_) newSet.insert(Cell(v.get_element(), valueMap.at(v.get_row_index())));
	column_.swap(newSet);
}

template<class Field_element_type>
inline void Set_column<Field_element_type>::add(Set_column &column)
{
	for (const Cell& v : column.column_){
		auto c = column_.find(v);
		if (c != column_.end()){
			c->get_element() += v.get_element();
			if (c->get_element() == 0) column_.erase(c);
		} else
			column_.insert(v);
	}
}

template<class Field_element_type>
inline Set_column<Field_element_type> &Set_column<Field_element_type>::operator=(Set_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	return *this;
}

template<class Friend_field_element_type>
inline void swap(Set_column<Friend_field_element_type>& col1, Set_column<Friend_field_element_type>& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // SETCOLUMN_H
