/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef B_Z2_UNORDERED_SET_COLUMN_H
#define B_Z2_UNORDERED_SET_COLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>
#include <algorithm>

#include "../../utilities/utilities.h"
#include "../z2_unordered_set_column.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Cell_type, class Row_access_option>
class Z2_unordered_set_boundary_column : public Z2_unordered_set_column<Cell_type,Row_access_option>
{
private:
	using Base = Z2_unordered_set_column<Cell_type,Row_access_option>;
	using Base::operator+=;		//kinda ugly, so TODO: organize better

public:
	using Cell = typename Base::Cell;
	using Column_type = typename Base::Column_type;
	using iterator = typename Base::iterator;
	using const_iterator = typename Base::const_iterator;

	Z2_unordered_set_boundary_column();
	template<class Boundary_type>
	Z2_unordered_set_boundary_column(const Boundary_type& boundary);
	template<class Boundary_type>
	Z2_unordered_set_boundary_column(const Boundary_type& boundary, dimension_type dimension);
	template<class Row_container_type>
	Z2_unordered_set_boundary_column(index columnIndex, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Z2_unordered_set_boundary_column(index columnIndex, const Boundary_type& boundary, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Z2_unordered_set_boundary_column(index columnIndex, const Boundary_type& boundary, dimension_type dimension, Row_container_type &rowContainer);
	Z2_unordered_set_boundary_column(const Z2_unordered_set_boundary_column& column);
	Z2_unordered_set_boundary_column(const Z2_unordered_set_boundary_column& column, index columnIndex);
	template<class Row_container_type>
	Z2_unordered_set_boundary_column(const Z2_unordered_set_boundary_column& column, index columnIndex, Row_container_type &rowContainer);
	Z2_unordered_set_boundary_column(Z2_unordered_set_boundary_column&& column) noexcept;

	int get_pivot();
	void clear();
	void clear(index rowIndex);
	template<class Map_type>
	void reorder(Map_type& valueMap);

	Z2_unordered_set_boundary_column& operator+=(Z2_unordered_set_boundary_column const &column);
	friend Z2_unordered_set_boundary_column operator+(Z2_unordered_set_boundary_column column1,
											 Z2_unordered_set_boundary_column const& column2){
		column1 += column2;
		return column1;
	}

	Z2_unordered_set_boundary_column& operator*=(unsigned int v);
	friend Z2_unordered_set_boundary_column operator*(Z2_unordered_set_boundary_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Z2_unordered_set_boundary_column operator*(unsigned int const& v, Z2_unordered_set_boundary_column column){
		column *= v;
		return column;
	}

	Z2_unordered_set_boundary_column& operator=(Z2_unordered_set_boundary_column other);

	friend void swap(Z2_unordered_set_boundary_column& col1,
					 Z2_unordered_set_boundary_column& col2){
		swap(static_cast<Z2_unordered_set_column<Cell_type,Row_access_option>&>(col1),
			 static_cast<Z2_unordered_set_column<Cell_type,Row_access_option>&>(col2));
		std::swap(col1.pivotChanged_, col2.pivotChanged_);
		std::swap(col1.pivot_, col2.pivot_);
	}

private:
	bool pivotChanged_;
	int pivot_;

	void _insert_cell(index rowIndex, Column_type& column);
};

template<class Cell_type, class Row_access_option>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column()
	: Base(), pivotChanged_(false), pivot_(-1)
{}

template<class Cell_type, class Row_access_option>
template<class Boundary_type>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column(const Boundary_type &boundary)
	: Base(boundary),
	  pivotChanged_(false),
	  pivot_(boundary.begin() == boundary.end() ? -1 : *std::prev(boundary.end()))
{}

template<class Cell_type, class Row_access_option>
template<class Boundary_type>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column(const Boundary_type &boundary, dimension_type dimension)
	: Base(boundary, dimension),
	  pivotChanged_(false),
	  pivot_(boundary.begin() == boundary.end() ? -1 : *std::prev(boundary.end()))
{}

template<class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column(
		index columnIndex, Row_container_type &rowContainer)
	: Base(columnIndex, rowContainer), pivotChanged_(false), pivot_(-1)
{}

template<class Cell_type, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column(
		index columnIndex, const Boundary_type &boundary, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, rowContainer),
	  pivotChanged_(false),
	  pivot_(boundary.begin() == boundary.end() ? -1 : *std::prev(boundary.end()))
{}

template<class Cell_type, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column(
		index columnIndex, const Boundary_type &boundary, dimension_type dimension, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, dimension, rowContainer),
	  pivotChanged_(false),
	  pivot_(boundary.begin() == boundary.end() ? -1 : *std::prev(boundary.end()))
{}

template<class Cell_type, class Row_access_option>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column(
		const Z2_unordered_set_boundary_column &column)
	: Base(static_cast<const Base&>(column)),
	  pivotChanged_(column.pivotChanged_),
	  pivot_(column.pivot_)
{}

template<class Cell_type, class Row_access_option>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column(
		const Z2_unordered_set_boundary_column &column, index columnIndex)
	: Base(static_cast<const Base&>(column), columnIndex),
	  pivotChanged_(column.pivotChanged_),
	  pivot_(column.pivot_)
{}

template<class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column(
		const Z2_unordered_set_boundary_column& column, index columnIndex, Row_container_type &rowContainer)
	: Base(static_cast<const Base&>(column), columnIndex, rowContainer),
	  pivotChanged_(column.pivotChanged_),
	  pivot_(column.pivot_)
{}

template<class Cell_type, class Row_access_option>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::Z2_unordered_set_boundary_column(
		Z2_unordered_set_boundary_column &&column) noexcept
	: Base(std::move(static_cast<Base&&>(column))),
	  pivotChanged_(std::exchange(column.pivotChanged_, 0)),
	  pivot_(std::exchange(column.pivot_, 0))
{}

template<class Cell_type, class Row_access_option>
inline int Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::get_pivot()
{
	if (pivotChanged_ && Base::column_.size() == 0){
		pivot_ = -1;
		pivotChanged_ = false;
	} else if (pivotChanged_) {
		pivot_ = 0;
		for (const Cell& c : Base::column_){
			if (static_cast<int>(c.get_row_index()) > pivot_)
				pivot_ = c.get_row_index();
		}
		pivotChanged_ = false;
	}

	return pivot_;
}

template<class Cell_type, class Row_access_option>
inline void Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::clear()
{
	Base::clear();
	pivot_ = -1;
	pivotChanged_ = false;
}

template<class Cell_type, class Row_access_option>
inline void Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::clear(index rowIndex)
{
	iterator it;
	if constexpr (Row_access_option::isActive_){
		it = Base::column_.find(Cell(Row_access_option::columnIndex_, rowIndex));
	} else {
		it = Base::column_.find(Cell(rowIndex));
	}
	if (it != Base::column_.end())
		Base::_delete_cell(it);
	if (static_cast<int>(rowIndex) == pivot_) pivotChanged_ = true;
}

template<class Cell_type, class Row_access_option>
template<class Map_type>
inline void Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::reorder(Map_type &valueMap)
{
	Base::reorder(valueMap);
	pivotChanged_ = true;
}

template<class Cell_type, class Row_access_option>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option> &
Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::operator+=(Z2_unordered_set_boundary_column const &column)
{
	for (const Cell& v : column.column_){
		int id = v.get_row_index();
		auto it = Base::column_.find(v);
		if (it != Base::column_.end()){
			Base::_delete_cell(it);
			if (id == pivot_) pivotChanged_ = true;
		} else {
			Base::_insert_cell(id);
			if (id > pivot_){
				pivot_ = id;
				pivotChanged_ = false;
			}
		}
	}

	return *this;
}

template<class Cell_type, class Row_access_option>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option> &
Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::operator*=(unsigned int v)
{
	if (v % 2 == 0)
		clear();

	return *this;
}

template<class Cell_type, class Row_access_option>
inline Z2_unordered_set_boundary_column<Cell_type,Row_access_option> &
Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::operator=(Z2_unordered_set_boundary_column other)
{
	Base::operator=(static_cast<Base&>(other));
	std::swap(pivotChanged_, other.pivotChanged_);
	std::swap(pivot_, other.pivot_);
	return *this;
}

template<class Cell_type, class Row_access_option>
inline void Z2_unordered_set_boundary_column<Cell_type,Row_access_option>::_insert_cell(
		index rowIndex, Column_type &column)
{
	if constexpr (Row_access_option::isActive_){
		auto it = column.try_emplace(Row_access_option::columnIndex_, rowIndex);
		Row_access_option::insert_cell(rowIndex, *it.first);
	} else {
		column.try_emplace(rowIndex);
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // B_Z2_UNORDERED_SET_COLUMN_H
