/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef B_Z2_VECTOR_COLUMN_H
#define B_Z2_VECTOR_COLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>

#include "../../utilities/utilities.h"
#include "../z2_vector_column.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Cell_type, class Row_access_option>
class Z2_vector_boundary_column : public Z2_vector_column<Cell_type,Row_access_option>
{
private:
	using Base = Z2_vector_column<Cell_type,Row_access_option>;
	using Base::operator+=;		//kinda ugly, so TODO: organize better

public:
	using Cell = typename Base::Cell;
	using Column_type = typename Base::Column_type;
	using iterator = typename Base::iterator;
	using const_iterator = typename Base::const_iterator;

	Z2_vector_boundary_column();
	template<class Boundary_type>
	Z2_vector_boundary_column(const Boundary_type& boundary);
	template<class Boundary_type>
	Z2_vector_boundary_column(const Boundary_type& boundary, dimension_type dimension);
	template<class Row_container_type>
	Z2_vector_boundary_column(index columnIndex, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Z2_vector_boundary_column(index columnIndex, const Boundary_type& boundary, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Z2_vector_boundary_column(index columnIndex, const Boundary_type& boundary, dimension_type dimension, Row_container_type &rowContainer);
	Z2_vector_boundary_column(const Z2_vector_boundary_column& column);
	Z2_vector_boundary_column(const Z2_vector_boundary_column& column, index columnIndex);
	Z2_vector_boundary_column(Z2_vector_boundary_column&& column) noexcept;

	std::vector<bool> get_content(int columnLength = -1);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	int get_pivot();
	void clear();
	void clear(index rowIndex);
	template<class Map_type>
	void reorder(Map_type& valueMap);

	Z2_vector_boundary_column& operator+=(Z2_vector_boundary_column &column);
	friend Z2_vector_boundary_column operator+(Z2_vector_boundary_column column1, Z2_vector_boundary_column& column2){
		column1 += column2;
		return column1;
	}

	Z2_vector_boundary_column& operator*=(unsigned int v);
	friend Z2_vector_boundary_column operator*(Z2_vector_boundary_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Z2_vector_boundary_column operator*(unsigned int const& v, Z2_vector_boundary_column column){
		column *= v;
		return column;
	}

	Z2_vector_boundary_column& operator=(Z2_vector_boundary_column other);

	friend void swap(Z2_vector_boundary_column& col1, Z2_vector_boundary_column& col2){
		swap(static_cast<Z2_vector_column<Cell_type,Row_access_option>&>(col1),
			 static_cast<Z2_vector_column<Cell_type,Row_access_option>&>(col2));
		col1.erasedValues_.swap(col2.erasedValues_);
	}

private:
	std::unordered_set<unsigned int> erasedValues_;

	void _cleanValues();
};

template<class Cell_type, class Row_access_option>
inline Z2_vector_boundary_column<Cell_type,Row_access_option>::Z2_vector_boundary_column() : Base()
{}

template<class Cell_type, class Row_access_option>
template<class Boundary_type>
inline Z2_vector_boundary_column<Cell_type,Row_access_option>::Z2_vector_boundary_column(const Boundary_type &boundary)
	: Base(boundary)
{}

template<class Cell_type, class Row_access_option>
template<class Boundary_type>
inline Z2_vector_boundary_column<Cell_type,Row_access_option>::Z2_vector_boundary_column(const Boundary_type &boundary, dimension_type dimension)
	: Base(boundary, dimension)
{}

template<class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Z2_vector_boundary_column<Cell_type,Row_access_option>::Z2_vector_boundary_column(
		index columnIndex, Row_container_type &rowContainer)
	: Base(columnIndex, rowContainer)
{}

template<class Cell_type, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Z2_vector_boundary_column<Cell_type,Row_access_option>::Z2_vector_boundary_column(
		index columnIndex, const Boundary_type &boundary, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, rowContainer)
{}

template<class Cell_type, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Z2_vector_boundary_column<Cell_type,Row_access_option>::Z2_vector_boundary_column(
		index columnIndex, const Boundary_type &boundary, dimension_type dimension, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, dimension, rowContainer)
{}

template<class Cell_type, class Row_access_option>
inline Z2_vector_boundary_column<Cell_type,Row_access_option>::Z2_vector_boundary_column(
		const Z2_vector_boundary_column &column)
	: Base(static_cast<const Base&>(column)),
	  erasedValues_(column.erasedValues_)
{}

template<class Cell_type, class Row_access_option>
inline Z2_vector_boundary_column<Cell_type,Row_access_option>::Z2_vector_boundary_column(
		const Z2_vector_boundary_column &column, index columnIndex)
	: Base(static_cast<const Base&>(column), columnIndex),
	  erasedValues_(column.erasedValues_)
{}

template<class Cell_type, class Row_access_option>
inline Z2_vector_boundary_column<Cell_type,Row_access_option>::Z2_vector_boundary_column(Z2_vector_boundary_column &&column) noexcept
	: Base(std::move(static_cast<Base&&>(column))),
	  erasedValues_(std::move(column.erasedValues_))
{}

template<class Cell_type, class Row_access_option>
inline std::vector<bool> Z2_vector_boundary_column<Cell_type,Row_access_option>::get_content(int columnLength)
{
	_cleanValues();
	return Base::get_content(columnLength);
}

template<class Cell_type, class Row_access_option>
inline bool Z2_vector_boundary_column<Cell_type,Row_access_option>::is_non_zero(index rowIndex) const
{
	if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

	return Base::is_non_zero(rowIndex);
}

template<class Cell_type, class Row_access_option>
inline bool Z2_vector_boundary_column<Cell_type,Row_access_option>::is_empty()
{
	_cleanValues();
	return Base::column_.empty();
}

template<class Cell_type, class Row_access_option>
inline int Z2_vector_boundary_column<Cell_type,Row_access_option>::get_pivot()
{
	if (Base::column_.empty()) return -1;

	auto it = erasedValues_.find(Base::column_.back()->get_row_index());
	while (!Base::column_.empty() && it != erasedValues_.end()) {
		erasedValues_.erase(it);
		Base::_delete_cell(Base::column_.back());
		Base::column_.pop_back();
		if (!Base::column_.empty()) it = erasedValues_.find(Base::column_.back()->get_row_index());
	}

	if (Base::column_.empty()) return -1;

	return Base::column_.back()->get_row_index();
}

template<class Cell_type, class Row_access_option>
inline void Z2_vector_boundary_column<Cell_type,Row_access_option>::clear()
{
	Base::clear();
	erasedValues_.clear();
}

template<class Cell_type, class Row_access_option>
inline void Z2_vector_boundary_column<Cell_type,Row_access_option>::clear(index rowIndex)
{
	erasedValues_.insert(rowIndex);
}

template<class Cell_type, class Row_access_option>
template<class Map_type>
inline void Z2_vector_boundary_column<Cell_type,Row_access_option>::reorder(Map_type &valueMap)
{
	Column_type newColumn;
	for (Cell* v : Base::column_) {
		if (erasedValues_.find(v->get_row_index()) == erasedValues_.end()){
			v->set_row_index(valueMap[v->get_row_index()]);
			newColumn.push_back(v);
			if constexpr (Row_access_option::isActive_){
				Row_access_option::unlink(v);
			}
		} else {
			Base::_delete_cell(v);
		}
	}
	//all cells have to be deleted first, to avoid problem with insertion when row is a set
	if constexpr (Row_access_option::isActive_){
		for (Cell* cell : Base::column_) {
			Row_access_option::insert_cell(cell->get_row_index(), cell);
		}
	}
	std::sort(newColumn.begin(), newColumn.end(), [](const Cell* c1, const Cell* c2){return *c1 < *c2;});
	erasedValues_.clear();
	Base::column_.swap(newColumn);
}

template<class Cell_type, class Row_access_option>
inline Z2_vector_boundary_column<Cell_type,Row_access_option> &Z2_vector_boundary_column<Cell_type,Row_access_option>::operator+=(Z2_vector_boundary_column &column)
{
	_cleanValues();
	column._cleanValues();
	Base::operator+=(column);

	return *this;
}

template<class Cell_type, class Row_access_option>
inline Z2_vector_boundary_column<Cell_type,Row_access_option> &Z2_vector_boundary_column<Cell_type,Row_access_option>::operator*=(unsigned int v)
{
	Base::operator*=(v);

	if (Base::column_.empty())
		erasedValues_.clear();

	return *this;
}

template<class Cell_type, class Row_access_option>
inline Z2_vector_boundary_column<Cell_type,Row_access_option> &Z2_vector_boundary_column<Cell_type,Row_access_option>::operator=(Z2_vector_boundary_column other)
{
	Base::operator=(static_cast<Base&>(other));
	erasedValues_.swap(other.erasedValues_);
	return *this;
}

template<class Cell_type, class Row_access_option>
inline void Z2_vector_boundary_column<Cell_type,Row_access_option>::_cleanValues()
{
	if (erasedValues_.empty()) return;

	Column_type newColumn;
	for (Cell* v : Base::column_){
		if (erasedValues_.find(v->get_row_index()) == erasedValues_.end())
			newColumn.push_back(v);
		else
			Base::_delete_cell(v);
	}
	erasedValues_.clear();
	Base::column_.swap(newColumn);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // B_Z2_VECTOR_COLUMN_H
