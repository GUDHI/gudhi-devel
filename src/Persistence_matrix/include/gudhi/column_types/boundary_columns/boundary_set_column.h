/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef B_SET_COLUMN_H
#define B_SET_COLUMN_H

#include <iostream>
#include <list>
#include <set>

#include "../../utilities/utilities.h"
#include "../set_column.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Cell_type, class Row_access_option>
class Set_boundary_column : public Set_column<Field_element_type,Cell_type,Row_access_option>
{
private:
	using Base = Set_column<Field_element_type,Cell_type,Row_access_option>;
	using Base::operator+=;				//kinda ugly, so TODO: organize better
	using Base::multiply_and_add;		//kinda ugly, so TODO: organize better

public:
	using Cell = typename Base::Cell;
	using Column_type = typename Base::Column_type;
	using iterator = typename Base::iterator;
	using const_iterator = typename Base::const_iterator;

	Set_boundary_column();
	template<class Boundary_type>
	Set_boundary_column(const Boundary_type& boundary);
	template<class Boundary_type>
	Set_boundary_column(const Boundary_type& boundary, dimension_type dimension);
	template<class Row_container_type>
	Set_boundary_column(index columnIndex, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Set_boundary_column(index columnIndex, const Boundary_type& boundary, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Set_boundary_column(index columnIndex, const Boundary_type& boundary, dimension_type dimension, Row_container_type &rowContainer);
	Set_boundary_column(const Set_boundary_column& column);
	Set_boundary_column(const Set_boundary_column& column, index columnIndex);
	Set_boundary_column(Set_boundary_column&& column) noexcept;

	int get_pivot() const;
	Field_element_type get_pivot_value() const;
	using Base::clear;
	void clear(index rowIndex);

	Set_boundary_column& operator+=(Set_boundary_column const &column);
	friend Set_boundary_column operator+(Set_boundary_column column1, Set_boundary_column const& column2){
		column1 += column2;
		return column1;
	}
	friend Set_boundary_column operator*(Set_boundary_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Set_boundary_column operator*(unsigned int const& v, Set_boundary_column column){
		column *= v;
		return column;
	}

	//this = v * this + column
	Set_boundary_column& multiply_and_add(const Field_element_type& v, const Set_boundary_column& column);
	//this = this + column * v
	Set_boundary_column& multiply_and_add(const Set_boundary_column& column, const Field_element_type& v);

	Set_boundary_column& operator=(Set_boundary_column other);

	friend void swap(Set_boundary_column& col1, Set_boundary_column& col2){
		swap(static_cast<Set_column<Field_element_type,Cell_type,Row_access_option>&>(col1),
			 static_cast<Set_column<Field_element_type,Cell_type,Row_access_option>&>(col2));
	}

private:
	void _insert_cell(const Field_element_type& value, index rowIndex, Column_type& column);
};

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Set_boundary_column() : Base()
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Boundary_type>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Set_boundary_column(const Boundary_type &boundary)
	: Base(boundary)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Boundary_type>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Set_boundary_column(const Boundary_type &boundary, dimension_type dimension)
	: Base(boundary, dimension)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Set_boundary_column(
		index columnIndex, Row_container_type &rowContainer)
	: Base(columnIndex, rowContainer)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Set_boundary_column(
		index columnIndex, const Boundary_type& boundary, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, rowContainer)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Set_boundary_column(
		index columnIndex, const Boundary_type& boundary, dimension_type dimension, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, dimension, rowContainer)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Set_boundary_column(
		const Set_boundary_column &column)
	: Base(static_cast<const Base&>(column))
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Set_boundary_column(
		const Set_boundary_column& column, index columnIndex)
	: Base(static_cast<const Base&>(column), columnIndex)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::Set_boundary_column(
		Set_boundary_column &&column) noexcept
	: Base(std::move(static_cast<Base&&>(column)))
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline int Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::get_pivot() const
{
	if (Base::column_.empty()) return -1;
	return Base::column_.rbegin()->get_row_index();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Field_element_type Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::get_pivot_value() const
{
	if (Base::column_.empty()) return 0;
	return Base::column_.rbegin()->get_element();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::clear(index rowIndex)
{
	iterator it;
	if constexpr (Row_access_option::isActive_){
		it = Base::column_.find(Cell(0, Row_access_option::columnIndex_, rowIndex));
	} else {
		it = Base::column_.find(Cell(0, rowIndex));
	}
	if (it != Base::column_.end())
		Base::_delete_cell(it);
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option> &
Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::operator+=(Set_boundary_column const &column)
{
	Base::operator+=(column);

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option> &
Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::multiply_and_add(const Field_element_type& val, const Set_boundary_column& column)
{
	Base::multiply_and_add(val, column);

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option> &
Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::multiply_and_add(const Set_boundary_column& column, const Field_element_type& val)
{
	Base::multiply_and_add(column, val);

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_boundary_column<Field_element_type,Cell_type,Row_access_option> &Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::operator=(Set_boundary_column other)
{
	Base::operator=(static_cast<Base&>(other));
	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Set_boundary_column<Field_element_type,Cell_type,Row_access_option>::_insert_cell(
		const Field_element_type &value, index rowIndex, Column_type &column)
{
	if constexpr (Row_access_option::isActive_){
		auto it = column.emplace_hint(column.end(), value, Row_access_option::columnIndex_, rowIndex);
		Row_access_option::insert_cell(rowIndex, *it);
	} else {
		column.emplace_hint(column.end(), value, rowIndex);
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // B_SET_COLUMN_H
