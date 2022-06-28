/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef LISTCOLUMN_H
#define LISTCOLUMN_H

#include <iostream>
#include <list>
#include <vector>

#include "../utilities.h"
#include "../Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type = Zp_field_element<11> >
class List_column
{
public:
	using Cell = Base_cell<Field_element_type>;

	List_column();
	List_column(std::vector<index>& rowIndices, std::vector<unsigned int>& values);
	List_column(List_column& column);
	List_column(List_column&& column) noexcept;

//	void get_content(boundary_type& container);
	bool isNonZero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);
	void add(List_column& column);

	List_column& operator=(List_column other);

	template<class Friend_field_element_type>
	friend void swap(List_column<Friend_field_element_type>& col1,
					 List_column<Friend_field_element_type>& col2);

private:
	int dim_;
	std::list<Cell> column_;
};

template<class Field_element_type>
inline List_column<Field_element_type>::List_column() : dim_(0)
{}

template<class Field_element_type>
inline List_column<Field_element_type>::List_column(std::vector<index>& rowIndices, std::vector<unsigned int>& values)
	: dim_(rowIndices.size() == 0 ? 0 : rowIndices.size() - 1)
{
	for (unsigned int i = 0; i < rowIndices.size(); i++){
		column_.push_back(Cell(values.at(i), rowIndices[i]));
	}
}

template<class Field_element_type>
inline List_column<Field_element_type>::List_column(List_column &column)
	: dim_(column.dim_),
	  column_(column.column_)
{}

template<class Field_element_type>
inline List_column<Field_element_type>::List_column(List_column &&column) noexcept
	: dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

//template<class Field_element_type>
//inline void List_column<Field_element_type>::get_content(boundary_type &container)
//{
//	std::copy(column_.begin(), column_.end(), std::back_inserter(container));
//}

template<class Field_element_type>
inline bool List_column<Field_element_type>::isNonZero(index rowIndex) const
{
	for (Cell v : column_){
		if (v.get_row_index() == rowIndex) return true;
	}
	return false;
}

template<class Field_element_type>
inline bool List_column<Field_element_type>::is_empty()
{
	return column_.empty();
}

template<class Field_element_type>
inline dimension_type List_column<Field_element_type>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type>
inline int List_column<Field_element_type>::get_pivot()
{
	if (column_.empty()) return -1;

	return column_.back();
}

template<class Field_element_type>
inline void List_column<Field_element_type>::clear()
{
	column_.clear();
}

template<class Field_element_type>
inline void List_column<Field_element_type>::clear(index rowIndex)
{
	auto it = column_.begin();
	while (it != column_.end() && it->get_row_index() != rowIndex) it++;
	if (it != column_.end()) column_.erase(it);
}

template<class Field_element_type>
inline void List_column<Field_element_type>::reorder(std::vector<index> &valueMap)
{
	typename std::list<Cell>::iterator it = column_.begin();
	while (it != column_.end()) {
		it->setRowIndex(valueMap.at(it->get_row_index()));
		it++;
	}
	column_.sort();
}

template<class Field_element_type>
inline void List_column<Field_element_type>::add(List_column &column)
{
	if (column.is_empty()) return;
	if (column_.empty()){
		std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
		return;
	}

	typename std::list<Cell>::iterator itToAdd = column.column_.begin();
	typename std::list<Cell>::iterator itTarget = column_.begin();
	index curRowToAdd = itToAdd->get_row_index();
	index curRowTarget = itTarget->get_row_index();

	while (itToAdd != column.column_.end() && itTarget != column_.end())
	{
		if (itToAdd != column.column_.end() && itTarget != column_.end()){
			if (curRowToAdd == curRowTarget){
				itTarget->get_element() =+ itToAdd->get_element();
				if (itTarget->get_element() == 0) column_.erase(itTarget++);
				else itTarget++;
				itToAdd++;
			} else if (curRowToAdd < curRowTarget){
				column_.insert(itTarget, Cell(itToAdd->get_element(), curRowToAdd));
				itToAdd++;
			} else {
				itTarget++;
			}
		}

		curRowToAdd = itToAdd->get_row_index();
		curRowTarget = itTarget->get_row_index();
	}

	while (itToAdd != column.column_.end()){
		curRowToAdd = itToAdd->get_row_index();
		if (itToAdd != column.column_.end()){
			column_.push_back(Cell(itToAdd->get_element(), curRowToAdd));
			itToAdd++;
		}
	}
}

template<class Field_element_type>
inline List_column<Field_element_type> &List_column<Field_element_type>::operator=(List_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	return *this;
}

template<class Friend_field_element_type>
inline void swap(List_column<Friend_field_element_type>& col1,
				 List_column<Friend_field_element_type>& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // LISTCOLUMN_H
