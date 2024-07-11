/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-24 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file row_access.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Row_access class and
 * @ref Gudhi::persistence_matrix::Dummy_row_access structure.
 */

#ifndef PM_ROW_ACCESS_H
#define PM_ROW_ACCESS_H

#include <utility>  //std::swap

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Row_access, if the row access is not enabled.
 */
struct Dummy_row_access 
{
  friend void swap([[maybe_unused]] Dummy_row_access& d1, [[maybe_unused]] Dummy_row_access& d2) {}

  Dummy_row_access() {}
  template <typename index, class Row_container_type>
  Dummy_row_access([[maybe_unused]] index columnIndex, [[maybe_unused]] Row_container_type& rows) {}
};

/**
 * @class Row_access row_access.h gudhi/Persistence_matrix/columns/row_access.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the row access for the columns.
 * 
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Row_access 
{
 public:
  using index = typename Master_matrix::index;                            /**< @ref MatIdx index type. */
  using id_index = typename Master_matrix::id_index;                      /**< @ref IDIdx index type. */
  using Cell_type = typename Master_matrix::Cell_type;                    /**< @ref Cell. */
  using Row_container_type = typename Master_matrix::row_container_type;  /**< Type of the row container. */

  /**
   * @brief Default constructor. Sets the column index to -1 and the row container to nullptr.
   * Should only be used by dummy columns.
   */
  Row_access();
  /**
   * @brief Constructor setting the column index and the row container by the given values.
   * 
   * @param columnIndex Column index to store.
   * @param rows Pointer to the row container.
   */
  Row_access(index columnIndex, Row_container_type* rows);
  /**
   * @brief Move constructor.
   * 
   * @param other Column to move.
   */
  Row_access(Row_access&& other) noexcept;

  /**
   * @brief Inserts the given cell at the given row index.
   * 
   * @param rowIndex @ref rowindex "Row index" of the cell.
   * @param cell Pointer to the cell to insert.
   */
  void insert_cell(id_index rowIndex, Cell_type* cell);
  /**
   * @brief Removes the given cell from its row.
   * 
   * @param cell Pointer to the cell to remove.
   */
  void unlink(Cell_type* cell);
  /**
   * @brief If @ref PersistenceMatrixOptions::has_intrusive_rows is false, updates the copy of the cell in its row. 
   * Otherwise does nothing.
   *
   * If the rows are intrusive, only a pointer of the cell is stored and therefore any update on the cell (value
   * or column index) is automatically forwarded. But for non intrusive rows, any update has to be pushed explicitly.
   * 
   * @param cell Cell to update.
   */
  void update_cell(const Cell_type& cell);
  /**
   * @brief Returns the @ref MatIdx column index.
   * 
   * @return The @ref MatIdx column index.
   */
  index get_column_index() const;

  /**
   * @brief Swap operator.
   */
  friend void swap(Row_access& r1, Row_access& r2) {
    std::swap(r1.rows_, r2.rows_);
    std::swap(r1.columnIndex_, r2.columnIndex_);
  }

  // void set_rows(Row_container_type *rows);

 protected:
  index columnIndex_;         /**< Column index. */
  Row_container_type* rows_;  /**< Row container. Be careful to not destroy before the columns. */

 private:
  using base_hook_matrix_row = typename Master_matrix::base_hook_matrix_row;
};

template <class Master_matrix>
inline Row_access<Master_matrix>::Row_access() : columnIndex_(-1), rows_(nullptr) 
{}

template <class Master_matrix>
inline Row_access<Master_matrix>::Row_access(index columnIndex, Row_container_type* rows)
    : columnIndex_(columnIndex), rows_(rows) 
{}

template <class Master_matrix>
inline Row_access<Master_matrix>::Row_access(Row_access&& other) noexcept
    : columnIndex_(std::exchange(other.columnIndex_, 0)), rows_(other.rows_) 
{}

template <class Master_matrix>
inline void Row_access<Master_matrix>::insert_cell(id_index rowIndex, Cell_type* cell) 
{
  if (rows_ == nullptr) return;

  if constexpr (!Master_matrix::Option_list::has_removable_rows) {
    if (rows_->size() < rowIndex + 1) rows_->resize(rowIndex + 1);
  }

  // if has_removable_rows should op[] create non existing entry? If not, use try_emplace()
  if constexpr (Master_matrix::Option_list::has_intrusive_rows) {
    rows_->operator[](rowIndex).push_back(*cell);
  } else {
    rows_->operator[](rowIndex).insert(*cell);
  }
}

template <class Master_matrix>
inline void Row_access<Master_matrix>::unlink(Cell_type* cell) 
{
  if (rows_ == nullptr) return;

  if constexpr (Master_matrix::Option_list::has_intrusive_rows) {
    cell->base_hook_matrix_row::unlink();
  } else {
    if constexpr (Master_matrix::Option_list::has_removable_rows) {
      auto it = rows_->find(cell->get_row_index());
      it->second.erase(*cell);
    } else {
      rows_->operator[](cell->get_row_index()).erase(*cell);
    }
  }
}

template <class Master_matrix>
inline void Row_access<Master_matrix>::update_cell(const Cell_type& cell) 
{
  if constexpr (!Master_matrix::Option_list::has_intrusive_rows) {
    if (rows_ == nullptr) return;
    auto& row = rows_->at(cell.get_row_index());
    auto it = row.find(cell);
    it = row.erase(it);
    row.insert(it, cell);
  }
}

template <class Master_matrix>
inline typename Row_access<Master_matrix>::index Row_access<Master_matrix>::get_column_index() const 
{
  return columnIndex_;
}

// template<class Master_matrix>
// inline void Row_access<Master_matrix>::set_rows(Row_container_type *rows)
// {
// 	rows_ = rows;
// }

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_ROW_ACCESS_H
