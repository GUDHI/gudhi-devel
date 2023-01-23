/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef B_UNORDERED_SET_COLUMN_H
#define B_UNORDERED_SET_COLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>
#include <algorithm>

#include "../../utilities/utilities.h"
#include "../unordered_set_column.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Cell_type, class Row_access_option>
class Unordered_set_boundary_column : public Unordered_set_column<Field_element_type,Cell_type,Row_access_option>
{
private:
	using Base = Unordered_set_column<Field_element_type,Cell_type,Row_access_option>;

public:
	using Cell = typename Base::Cell;
	using Column_type = typename Base::Column_type;
	using iterator = typename Base::iterator;
	using const_iterator = typename Base::const_iterator;

	Unordered_set_boundary_column();
	template<class Boundary_type>
	Unordered_set_boundary_column(const Boundary_type& boundary);
	template<class Boundary_type>
	Unordered_set_boundary_column(const Boundary_type& boundary, dimension_type dimension);
	template<class Row_container_type>
	Unordered_set_boundary_column(index columnIndex, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Unordered_set_boundary_column(index columnIndex, const Boundary_type& boundary, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Unordered_set_boundary_column(index columnIndex, const Boundary_type& boundary, dimension_type dimension, Row_container_type &rowContainer);
	Unordered_set_boundary_column(const Unordered_set_boundary_column& column);
	Unordered_set_boundary_column(const Unordered_set_boundary_column& column, index columnIndex);
	Unordered_set_boundary_column(Unordered_set_boundary_column&& column) noexcept;

	int get_pivot();
	Field_element_type get_pivot_value();
	void clear();
	void clear(index rowIndex);
	template<class Map_type>
	void reorder(Map_type& valueMap);

	Unordered_set_boundary_column& operator+=(Unordered_set_boundary_column const &column);
	friend Unordered_set_boundary_column operator+(Unordered_set_boundary_column column1, Unordered_set_boundary_column const& column2){
		column1 += column2;
		return column1;
	}
	Unordered_set_boundary_column& operator*=(unsigned int v);
	friend Unordered_set_boundary_column operator*(Unordered_set_boundary_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Unordered_set_boundary_column operator*(unsigned int const& v, Unordered_set_boundary_column column){
		column *= v;
		return column;
	}

	Unordered_set_boundary_column& operator=(Unordered_set_boundary_column other);

	friend void swap(Unordered_set_boundary_column& col1,
					 Unordered_set_boundary_column& col2){
		swap(static_cast<Unordered_set_column<Field_element_type,Cell_type,Row_access_option>&>(col1),
			 static_cast<Unordered_set_column<Field_element_type,Cell_type,Row_access_option>&>(col2));
		std::swap(col1.pivotChanged_, col2.pivotChanged_);
		std::swap(col1.pivot_, col2.pivot_);
	}

private:
	bool pivotChanged_;
	Cell pivot_;

	void _insert_cell(const Field_element_type& value, index rowIndex, Column_type& column);
};

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Unordered_set_boundary_column()
	: Base(), pivotChanged_(false)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Boundary_type>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Unordered_set_boundary_column(const Boundary_type &boundary)
	: Base(boundary),
	  pivotChanged_(false),
	  pivot_(boundary.empty() ? Cell() : Cell(boundary.rbegin()->second, boundary.rbegin()->first))
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Boundary_type>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Unordered_set_boundary_column(const Boundary_type &boundary, dimension_type dimension)
	: Base(boundary, dimension),
	  pivotChanged_(false),
	  pivot_(boundary.empty() ? Cell() : Cell(boundary.rbegin()->second, boundary.rbegin()->first))
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Unordered_set_boundary_column(
		index columnIndex, Row_container_type &rowContainer)
	: Base(columnIndex, rowContainer), pivotChanged_(false)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Unordered_set_boundary_column(
		index columnIndex, const Boundary_type& boundary, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, rowContainer),
	  pivotChanged_(false)
{
	if (!boundary.empty()){
		if constexpr (Row_access_option::isActive_){
			pivot_ = Cell(boundary.rbegin()->second, columnIndex, boundary.rbegin()->first);
		} else {
			pivot_ = Cell(boundary.rbegin()->second, boundary.rbegin()->first);
		}
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Unordered_set_boundary_column(
		index columnIndex, const Boundary_type& boundary, dimension_type dimension, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, dimension, rowContainer),
	  pivotChanged_(false)
{
	if (!boundary.empty()){
		if constexpr (Row_access_option::isActive_){
			pivot_ = Cell(boundary.rbegin()->second, columnIndex, boundary.rbegin()->first);
		} else {
			pivot_ = Cell(boundary.rbegin()->second, boundary.rbegin()->first);
		}
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Unordered_set_boundary_column(const Unordered_set_boundary_column &column)
	: Base(static_cast<const Base&>(column)),
	  pivotChanged_(column.pivotChanged_),
	  pivot_(column.pivot_)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Unordered_set_boundary_column(
		const Unordered_set_boundary_column& column, index columnIndex)
	: Base(static_cast<const Base&>(column), columnIndex),
	  pivotChanged_(column.pivotChanged_),
	  pivot_(column.pivot_)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Unordered_set_boundary_column(Unordered_set_boundary_column &&column) noexcept
	: Base(std::move(static_cast<Base&&>(column))),
	  pivotChanged_(std::exchange(column.pivotChanged_, 0)),
	  pivot_(std::move(column.pivot_))
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline int Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::get_pivot()
{
	if (pivotChanged_){
		pivot_ = Base::column_.size() == 0 ?
					Cell()
				  : *std::max_element(Base::column_.begin(), Base::column_.end());
		pivotChanged_ = false;
	}

	if (pivot_.get_element() == 0u) return -1;
	return pivot_.get_row_index();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Field_element_type Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::get_pivot_value()
{
	if (pivotChanged_){
		pivot_ = Base::column_.size() == 0 ?
					Cell()
				  : *std::max_element(Base::column_.begin(), Base::column_.end());
		pivotChanged_ = false;
	}

	return pivot_.get_element();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::clear()
{
	if constexpr (Row_access_option::isActive_){
		for (const Cell& cell : Base::column_)
			Row_access_option::unlink(cell);
	}
	Base::column_.clear();
	pivot_ = Cell();
	pivotChanged_ = false;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::clear(index rowIndex)
{
	iterator it;
	if constexpr (Row_access_option::isActive_){
		it = Base::column_.find(Cell(0, Row_access_option::columnIndex_, rowIndex));
	} else {
		it = Base::column_.find(Cell(0, rowIndex));
	}
	if (it != Base::column_.end())
		Base::_delete_cell(it);
	if (rowIndex == pivot_.get_row_index()) pivotChanged_ = true;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Map_type>
inline void Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::reorder(Map_type &valueMap)
{
	Base::reorder(valueMap);
	pivotChanged_ = true;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option> &
Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::operator+=(Unordered_set_boundary_column const &column)
{
	for (const Cell& v : column.column_){
		auto c = Base::column_.find(v);
		if (c != Base::column_.end()){
			index r = c->get_row_index();
			Field_element_type coef = c->get_element();
			coef += v.get_element();
			Base::_delete_cell(c);
			if (coef != 0u) Base::_insert_cell(coef, r);
			else if (v == pivot_) pivotChanged_ = true;
		} else {
			Base::_insert_cell(v.get_element(), v.get_row_index());
			if (pivot_ < v){
				pivot_ = v;
				pivotChanged_ = false;
			}
		}
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option> &
Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::operator*=(unsigned int v)
{
	Base::operator*=(v);

	if (Base::column_.empty()){
		pivot_ = Cell();
		pivotChanged_ = false;
	} else {
		pivot_.get_element() *= v;
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option> &
Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::operator=(Unordered_set_boundary_column other)
{
	Base::operator=(static_cast<Base&>(other));
	std::swap(pivotChanged_, other.pivotChanged_);
	std::swap(pivot_, other.pivot_);
	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Unordered_set_boundary_column<Field_element_type,Cell_type,Row_access_option>::_insert_cell(
		const Field_element_type &value, index rowIndex, Column_type &column)
{
	if constexpr (Row_access_option::isActive_){
		auto it = column.emplace(value, Row_access_option::columnIndex_, rowIndex);
		Row_access_option::insert_cell(rowIndex, *it.first);
	} else {
		column.emplace(value, rowIndex);
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // B_UNORDERED_SET_COLUMN_H
