/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ROW_ACCESS_H
#define ROW_ACCESS_H

#include <utility>
#include <cassert>

#include "../utilities/utilities.h"
#include "cell.h"					//to recognize the alias "base_hook_matrix_row"

namespace Gudhi {
namespace persistence_matrix {

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
class Row_access
{
public:
	Row_access();
	Row_access(index columnIndex, Row_container_type& rows);
	Row_access(Row_access&& other) noexcept;

	void insert_cell(index rowIndex, Cell_type *cell);
	void insert_cell(index rowIndex, const Cell_type &cell);
	void unlink(Cell_type *cell);
	void unlink(const Cell_type &cell);
	void update_cell(const Cell_type &cell);
	index get_column_index() const;

	friend void swap(Row_access& r1, Row_access& r2){
		std::swap(r1.rows_, r2.rows_);
		std::swap(r1.columnIndex_, r2.columnIndex_);
	}

	void set_rows(Row_container_type *rows);

protected:
	index columnIndex_;
	Row_container_type* rows_;					//be carefull to not destroy container before columns
	static constexpr bool isActive_ = true;
};

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::Row_access() : rows_(nullptr)
{}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::Row_access(index columnIndex, Row_container_type &rows)
	: columnIndex_(columnIndex), rows_(&rows)
{}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::Row_access(Row_access &&other) noexcept
	: columnIndex_(std::exchange(other.columnIndex_, 0)),
	  rows_(other.rows_)
{}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::insert_cell(index rowIndex, Cell_type *cell)
{
	if (rows_ == nullptr) return;

	if constexpr (!hasRemovableRows){
		if (rows_->size() < rowIndex + 1)
			rows_->resize(rowIndex + 1);
	}	//if hasRemovableRows op[] should create non existing entry? If not, use try_emplace()

	if constexpr (isIntrusive) {
		rows_->operator[](rowIndex).push_back(*cell);
	} else {
		rows_->operator[](rowIndex).insert(*cell);
	}
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::insert_cell(index rowIndex, const Cell_type &cell)
{
	static_assert(!isIntrusive, "Cannot insert const cell in intrusive container.");

	if (rows_ == nullptr) return;

	if constexpr (hasRemovableRows){
		auto res = rows_->try_emplace(rowIndex);
		res.first->second.insert(cell);
	} else {
		if (rows_->size() < rowIndex + 1)
			rows_->resize(rowIndex + 1);
		rows_->operator[](rowIndex).insert(cell);
	}
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::unlink(Cell_type *cell)
{
	if (rows_ == nullptr) return;

	if constexpr (isIntrusive){
		cell->base_hook_matrix_row::unlink();
	} else {
		if constexpr (hasRemovableRows){
			auto it = rows_->find(cell->get_row_index());
			it->second.erase(*cell);
		} else {
			rows_->operator[](cell->get_row_index()).erase(*cell);
		}
	}
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::unlink(const Cell_type &cell)
{
	static_assert(!isIntrusive, "Cannot unlink const cell from intrusive container.");

	if (rows_ == nullptr) return;

	if constexpr (hasRemovableRows){
		auto it = rows_->find(cell.get_row_index());
		it->second.erase(cell);
	} else {
		rows_->operator[](cell.get_row_index()).erase(cell);
	}
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::update_cell(const Cell_type &cell)
{
	if constexpr (!isIntrusive){
		if (rows_ == nullptr) return;
		auto& row = rows_->at(cell.get_row_index());
		auto it = row.find(cell);
		it = row.erase(it);
		row.insert(it, cell);
	}
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline index Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::get_column_index() const
{
	return columnIndex_;
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableRows>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableRows>::set_rows(Row_container_type *rows)
{
	rows_ = rows;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // ROW_ACCESS_H
