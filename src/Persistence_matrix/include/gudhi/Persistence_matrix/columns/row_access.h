/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_ROW_ACCESS_H
#define PM_ROW_ACCESS_H

#include <utility>	//std::swap

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_row_access{
	friend void swap([[maybe_unused]] Dummy_row_access& d1, [[maybe_unused]] Dummy_row_access& d2){}

	Dummy_row_access(){}
	template<typename index, class Row_container_type>
	Dummy_row_access([[maybe_unused]] index columnIndex, [[maybe_unused]] Row_container_type& rows){}
	Dummy_row_access([[maybe_unused]] Dummy_row_access&& other) noexcept{}
};

template<class Master_matrix>
class Row_access
{
public:
	using index = typename Master_matrix::index;
	using Cell_type = typename Master_matrix::Cell_type;
	using Row_container_type = typename Master_matrix::row_container_type;

	Row_access();
	Row_access(index columnIndex, Row_container_type& rows);
	Row_access(Row_access&& other) noexcept;

	void insert_cell(index rowIndex, Cell_type *cell);
	void insert_cell(index rowIndex, const Cell_type &cell);	//still used??
	void unlink(Cell_type *cell);
	void unlink(const Cell_type &cell);							//still used??
	void update_cell(const Cell_type &cell);
	index get_column_index() const;

	friend void swap(Row_access& r1, Row_access& r2){
		std::swap(r1.rows_, r2.rows_);
		std::swap(r1.columnIndex_, r2.columnIndex_);
	}

	void set_rows(Row_container_type *rows);

protected:
	index columnIndex_;
	Row_container_type* rows_;			//be carefull to not destroy container before columns

private:
	using base_hook_matrix_row = typename Master_matrix::base_hook_matrix_row;
};

template<class Master_matrix>
inline Row_access<Master_matrix>::Row_access() : columnIndex_(-1), rows_(nullptr)
{}

template<class Master_matrix>
inline Row_access<Master_matrix>::Row_access(index columnIndex, Row_container_type &rows)
	: columnIndex_(columnIndex), rows_(&rows)
{}

template<class Master_matrix>
inline Row_access<Master_matrix>::Row_access(Row_access &&other) noexcept
	: columnIndex_(std::exchange(other.columnIndex_, 0)),
	  rows_(other.rows_)
{}

template<class Master_matrix>
inline void Row_access<Master_matrix>::insert_cell(index rowIndex, Cell_type *cell)
{
	if (rows_ == nullptr) return;

	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		if (rows_->size() < rowIndex + 1)
			rows_->resize(rowIndex + 1);
	}
	
	//if has_removable_rows should op[] create non existing entry? If not, use try_emplace()
	if constexpr (Master_matrix::Option_list::has_intrusive_rows) {
		rows_->operator[](rowIndex).push_back(*cell);
	} else {
		rows_->operator[](rowIndex).insert(*cell);
	}
}

template<class Master_matrix>
inline void Row_access<Master_matrix>::insert_cell(index rowIndex, const Cell_type &cell)
{
	static_assert(!Master_matrix::Option_list::has_intrusive_rows, "Cannot insert const cell in intrusive container.");

	if (rows_ == nullptr) return;

	if constexpr (Master_matrix::Option_list::has_removable_rows){
		auto res = rows_->try_emplace(rowIndex);
		res.first->second.insert(cell);
	} else {
		if (rows_->size() < rowIndex + 1)
			rows_->resize(rowIndex + 1);
		rows_->operator[](rowIndex).insert(cell);
	}
}

template<class Master_matrix>
inline void Row_access<Master_matrix>::unlink(Cell_type *cell)
{
	if (rows_ == nullptr) return;

	if constexpr (Master_matrix::Option_list::has_intrusive_rows){
		cell->base_hook_matrix_row::unlink();
	} else {
		if constexpr (Master_matrix::Option_list::has_removable_rows){
			auto it = rows_->find(cell->get_row_index());
			it->second.erase(*cell);
		} else {
			rows_->operator[](cell->get_row_index()).erase(*cell);
		}
	}
}

template<class Master_matrix>
inline void Row_access<Master_matrix>::unlink(const Cell_type &cell)
{
	static_assert(!Master_matrix::Option_list::has_intrusive_rows, "Cannot unlink const cell from intrusive container.");

	if (rows_ == nullptr) return;

	if constexpr (Master_matrix::Option_list::has_removable_rows){
		auto it = rows_->find(cell.get_row_index());
		it->second.erase(cell);
	} else {
		rows_->operator[](cell.get_row_index()).erase(cell);
	}
}

template<class Master_matrix>
inline void Row_access<Master_matrix>::update_cell(const Cell_type &cell)
{
	if constexpr (!Master_matrix::Option_list::has_intrusive_rows){
		if (rows_ == nullptr) return;
		auto& row = rows_->at(cell.get_row_index());
		auto it = row.find(cell);
		it = row.erase(it);
		row.insert(it, cell);
	}
}

template<class Master_matrix>
inline typename Row_access<Master_matrix>::index Row_access<Master_matrix>::get_column_index() const
{
	return columnIndex_;
}

template<class Master_matrix>
inline void Row_access<Master_matrix>::set_rows(Row_container_type *rows)
{
	rows_ = rows;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_ROW_ACCESS_H
