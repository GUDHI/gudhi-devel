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
#include "cell.h"	//why??

namespace Gudhi {
namespace persistence_matrix {

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
class Row_access
{
public:
	Row_access(index columnIndex, Row_container_type& rows);
//	Row_access(const Row_access &toCopy);
	Row_access(Row_access&& other) noexcept;

	void insert_cell(index rowIndex, Cell_type *cell);
	void insert_cell(index rowIndex, const Cell_type &cell);
	void unlink(Cell_type *cell);
	void unlink(const Cell_type &cell);
	void update_cell(const Cell_type &cell);
	index get_column_index() const;

//	Row_access& operator=(Row_access other);
	friend void swap(Row_access& r1, Row_access& r2){
		std::swap(r1.rows_, r2.rows_);
		std::swap(r1.columnIndex_, r2.columnIndex_);
	}

protected:
	index columnIndex_;
	Row_container_type* rows_;					//be carefull to not destroy container before columns
	static constexpr bool isActive_ = true;

private:
	Row_access();	//for more explicit error message when compiling
};

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
inline Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableColumns>::Row_access()
{}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
inline Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableColumns>::Row_access(index columnIndex, Row_container_type &rows)
	: columnIndex_(columnIndex), rows_(&rows)
{}

//template<class Master_matrix>
//inline Row_access<Master_matrix>::Row_access(const Row_access &toCopy)
//	: columnIndex_(toCopy.columnIndex_), rows_(toCopy.rows_)
//{}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
inline Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableColumns>::Row_access(Row_access &&other) noexcept
	: columnIndex_(std::exchange(other.columnIndex_, 0)), rows_(std::move(other.rows_))
{}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableColumns>::insert_cell(index rowIndex, Cell_type *cell)
{
	if constexpr (hasRemovableColumns)
		rows_->try_emplace(rowIndex);
	else {
		if (rows_->size() < rowIndex + 1)
			rows_->resize(rowIndex + 1);
	}
	if constexpr (isIntrusive){
		rows_->at(rowIndex).push_back(*cell);
	} else {
		rows_->at(rowIndex).insert(*cell);
	}
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableColumns>::insert_cell(index rowIndex, const Cell_type &cell)
{
	static_assert(!isIntrusive, "Cannot insert const cell in intrusive container.");

	if constexpr (hasRemovableColumns)
		rows_->try_emplace(rowIndex);
	else {
		if (rows_->size() < rowIndex + 1)
			rows_->resize(rowIndex + 1);
	}
	rows_->at(rowIndex).insert(cell);
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableColumns>::unlink(Cell_type *cell)
{
	if constexpr (isIntrusive){
		cell->base_hook_matrix_row::unlink();
	} else {
		rows_->at(cell->get_row_index()).erase(*cell);
	}
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableColumns>::unlink(const Cell_type &cell)
{
	static_assert(!isIntrusive, "Cannot unlink const cell from intrusive container.");

	rows_->at(cell.get_row_index()).erase(cell);
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
inline void Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableColumns>::update_cell(const Cell_type &cell)
{
	if constexpr (!isIntrusive){
		rows_->at(cell.get_row_index()).erase(cell);
		rows_->at(cell.get_row_index()).insert(cell);
	}
}

template<class Row_container_type, class Cell_type, bool isIntrusive, bool hasRemovableColumns>
inline index Row_access<Row_container_type,Cell_type,isIntrusive,hasRemovableColumns>::get_column_index() const
{
	return columnIndex_;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // ROW_ACCESS_H
