/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef B_LIST_COLUMN_H
#define B_LIST_COLUMN_H

#include <iostream>
#include <list>
#include <vector>

#include "../../utilities/utilities.h"
#include "../list_column.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
class List_boundary_column : public List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>
{
private:
	using Base = List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>;

public:
	using Cell = typename Base::Cell;
	using Column_type = typename Base::Column_type;
	using iterator = typename Base::iterator;
	using const_iterator = typename Base::const_iterator;

	List_boundary_column();
	template<class Boundary_type>
	List_boundary_column(const Boundary_type& boundary);
	template<class Boundary_type>
	List_boundary_column(const Boundary_type& boundary, dimension_type dimension);
	template<class Row_container_type>
	List_boundary_column(index columnIndex, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	List_boundary_column(index columnIndex, const Boundary_type& boundary, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	List_boundary_column(index columnIndex, const Boundary_type& boundary, dimension_type dimension, Row_container_type &rowContainer);
	List_boundary_column(const List_boundary_column& column);
	List_boundary_column(const List_boundary_column& column, index columnIndex);
	List_boundary_column(List_boundary_column&& column) noexcept;

	int get_pivot() const;
	Field_element_type get_pivot_value() const;
	void clear();
	void clear(index rowIndex);
	template<class Map_type>
	void reorder(Map_type& valueMap);

	friend List_boundary_column operator+(List_boundary_column column1, List_boundary_column const& column2){
		column1 += column2;
		return column1;
	}
	friend List_boundary_column operator*(List_boundary_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend List_boundary_column operator*(unsigned int const& v, List_boundary_column column){
		column *= v;
		return column;
	}

	List_boundary_column& operator=(List_boundary_column other);

	friend void swap(List_boundary_column& col1, List_boundary_column& col2){
		swap(static_cast<List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>&>(col1),
			 static_cast<List_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>&>(col2));
	}
};

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_boundary_column() : Base()
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Boundary_type>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_boundary_column(const Boundary_type &boundary)
	: Base(boundary)
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Boundary_type>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_boundary_column(const Boundary_type &boundary, dimension_type dimension)
	: Base(boundary, dimension)
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Row_container_type>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_boundary_column(
		index columnIndex, Row_container_type &rowContainer)
	: Base(columnIndex, rowContainer)
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_boundary_column(
		index columnIndex, const Boundary_type& boundary, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, rowContainer)
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_boundary_column(
		index columnIndex, const Boundary_type& boundary, dimension_type dimension, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, dimension, rowContainer)
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_boundary_column(
		const List_boundary_column &column)
	: Base(static_cast<const Base&>(column))
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_boundary_column(
		const List_boundary_column& column, index columnIndex)
	: Base(static_cast<const Base&>(column), columnIndex)
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::List_boundary_column(
		List_boundary_column &&column) noexcept
	: Base(std::move(static_cast<Base&&>(column)))
{}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline int List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::get_pivot() const
{
	if (Base::column_.empty()) return -1;

	return Base::column_.back().get_row_index();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline Field_element_type List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::get_pivot_value() const
{
	if (Base::column_.empty()) return 0;

	return Base::column_.back().get_element();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline void List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::clear()
{
	if constexpr (Row_access_option::isActive_){
		for (Cell& cell : Base::column_)
			Row_access_option::unlink(&cell);
	}
	Base::column_.clear();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline void List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::clear(index rowIndex)
{
	auto it = Base::column_.begin();
	while (it != Base::column_.end() && it->get_row_index() != rowIndex) it++;
	if (it != Base::column_.end()) Base::_delete_cell(it);
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Map_type>
inline void List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::reorder(Map_type &valueMap)
{
	iterator it = Base::column_.begin();
	while (it != Base::column_.end()) {
		Cell* cell = &(*it);
		if constexpr (Row_access_option::isActive_) Row_access_option::unlink(cell);
		cell->set_row_index(valueMap[cell->get_row_index()]);
		if constexpr (Row_access_option::isActive_) Row_access_option::insert_cell(cell->get_row_index(), cell);
		it++;
	}
	Base::column_.sort();
}

template<class Field_element_type, class Cell_type, class Column_pairing_option, class Row_access_option>
inline List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option> &
List_boundary_column<Field_element_type,Cell_type,Column_pairing_option,Row_access_option>::operator=(
		List_boundary_column other)
{
	Base::operator=(static_cast<Base&>(other));
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // B_LIST_COLUMN_H
