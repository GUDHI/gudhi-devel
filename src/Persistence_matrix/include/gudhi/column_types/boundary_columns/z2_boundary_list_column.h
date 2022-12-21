/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef B_Z2_LIST_COLUMN_H
#define B_Z2_LIST_COLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>

#include "../../utilities/utilities.h"
#include "../z2_list_column.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Cell_type, class Column_pairing_option, class Row_access_option>
class Z2_list_boundary_column : public Z2_list_column<Cell_type,Column_pairing_option,Row_access_option>
{
private:
	using Base = Z2_list_column<Cell_type,Column_pairing_option,Row_access_option>;

public:
	using Cell = typename Base::Cell;
	using Column_type = typename Base::Column_type;
	using iterator = typename Base::iterator;
	using const_iterator = typename Base::const_iterator;

	Z2_list_boundary_column();
	template<class Boundary_type>
	Z2_list_boundary_column(const Boundary_type& boundary);
	template<class Boundary_type>
	Z2_list_boundary_column(const Boundary_type& boundary, dimension_type dimension);
	template<class Row_container_type>
	Z2_list_boundary_column(index columnIndex, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Z2_list_boundary_column(index columnIndex, const Boundary_type& boundary, Row_container_type &rowContainer);
	template<class Boundary_type, class Row_container_type>
	Z2_list_boundary_column(index columnIndex, const Boundary_type& boundary, dimension_type dimension, Row_container_type &rowContainer);
	Z2_list_boundary_column(const Z2_list_boundary_column& column);
	Z2_list_boundary_column(const Z2_list_boundary_column& column, index columnIndex);
	Z2_list_boundary_column(Z2_list_boundary_column&& column) noexcept;

	int get_pivot() const;
	void clear();
	void clear(index rowIndex);
	template<class Map_type>
	void reorder(Map_type& valueMap);

	friend Z2_list_boundary_column operator+(Z2_list_boundary_column column1, Z2_list_boundary_column const& column2){
		column1 += column2;
		return column1;
	}
	friend Z2_list_boundary_column operator*(Z2_list_boundary_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Z2_list_boundary_column operator*(unsigned int const& v, Z2_list_boundary_column column){
		column *= v;
		return column;
	}

	Z2_list_boundary_column& operator=(Z2_list_boundary_column other);

	friend void swap(Z2_list_boundary_column& col1, Z2_list_boundary_column& col2){
		swap(static_cast<Z2_list_column<Cell_type,Column_pairing_option,Row_access_option>&>(col1),
			 static_cast<Z2_list_column<Cell_type,Column_pairing_option,Row_access_option>&>(col2));
	}
};

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_list_boundary_column() : Base()
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Boundary_type>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_list_boundary_column(const Boundary_type &boundary)
	: Base(boundary)
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Boundary_type>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_list_boundary_column(const Boundary_type &boundary, dimension_type dimension)
	: Base(boundary, dimension)
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Row_container_type>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_list_boundary_column(
		index columnIndex, Row_container_type &rowContainer)
	: Base(columnIndex, rowContainer)
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_list_boundary_column(
		index columnIndex, const Boundary_type &boundary, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, rowContainer)
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Boundary_type, class Row_container_type>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_list_boundary_column(
		index columnIndex, const Boundary_type &boundary, dimension_type dimension, Row_container_type &rowContainer)
	: Base(columnIndex, boundary, dimension, rowContainer)
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_list_boundary_column(
		const Z2_list_boundary_column &column)
	: Base(static_cast<const Base&>(column))
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_list_boundary_column(
		const Z2_list_boundary_column &column, index columnIndex)
	: Base(static_cast<const Base&>(column), columnIndex)
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_list_boundary_column(
		Z2_list_boundary_column &&column) noexcept
	: Base(std::move(static_cast<Base&&>(column)))
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline int Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::get_pivot() const
{
	if (Base::column_.empty()) return -1;

	return Base::column_.back().get_row_index();
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline void Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::clear()
{
	if constexpr (Row_access_option::isActive_){
		for (Cell& cell : Base::column_)
			Row_access_option::unlink(&cell);
	}
	Base::column_.clear();
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline void Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::clear(index rowIndex)
{
	auto it = Base::column_.begin();
	while (it != Base::column_.end() && it->get_row_index() != rowIndex) it++;
	if (it != Base::column_.end()) Base::_delete_cell(it);
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Map_type>
inline void Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::reorder(Map_type &valueMap)
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

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option> &Z2_list_boundary_column<Cell_type,Column_pairing_option,Row_access_option>::operator=(Z2_list_boundary_column other)
{
	Base::operator=(static_cast<Base&>(other));
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // B_Z2_LIST_COLUMN_H
