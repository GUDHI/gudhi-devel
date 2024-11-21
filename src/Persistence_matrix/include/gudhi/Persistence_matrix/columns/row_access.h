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
  template <typename Index, class Row_container>
  Dummy_row_access([[maybe_unused]] Index columnIndex, [[maybe_unused]] Row_container& rows) {}
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
  using Index = typename Master_matrix::Index;                  /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;            /**< @ref IDIdx index type. */
  using Matrix_entry = typename Master_matrix::Matrix_entry;    /**< @ref Entry. */
  using Row_container = typename Master_matrix::Row_container;  /**< Type of the row container. */

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
  Row_access(Index columnIndex, Row_container* rows);
  /**
   * @brief Move constructor.
   * 
   * @param other Column to move.
   */
  Row_access(Row_access&& other) noexcept;

  /**
   * @brief Inserts the given entry at the given row index.
   * 
   * @param rowIndex @ref rowindex "Row index" of the entry.
   * @param entry Pointer to the entry to insert.
   */
  void insert_entry(ID_index rowIndex, Matrix_entry* entry);
  /**
   * @brief Removes the given entry from its row.
   * 
   * @param entry Pointer to the entry to remove.
   */
  void unlink(Matrix_entry* entry);
  /**
   * @brief If @ref PersistenceMatrixOptions::has_intrusive_rows is false, updates the copy of the entry in its row. 
   * Otherwise does nothing.
   *
   * If the rows are intrusive, only a pointer of the entry is stored and therefore any update on the entry (value
   * or column index) is automatically forwarded. But for non intrusive rows, any update has to be pushed explicitly.
   * 
   * @param entry Entry to update.
   */
  void update_entry(const Matrix_entry& entry);
  /**
   * @brief Returns the @ref MatIdx column index.
   * 
   * @return The @ref MatIdx column index.
   */
  Index get_column_index() const;

  /**
   * @brief Swap operator.
   */
  friend void swap(Row_access& r1, Row_access& r2) {
    std::swap(r1.rows_, r2.rows_);
    std::swap(r1.columnIndex_, r2.columnIndex_);
  }

 protected:
  Index columnIndex_;         /**< Column index. */
  Row_container* rows_;  /**< Row container. Be careful to not destroy before the columns. */

 private:
  using Base_hook_matrix_row = typename Master_matrix::Base_hook_matrix_row;
};

template <class Master_matrix>
inline Row_access<Master_matrix>::Row_access() : columnIndex_(-1), rows_(nullptr) 
{}

template <class Master_matrix>
inline Row_access<Master_matrix>::Row_access(Index columnIndex, Row_container* rows)
    : columnIndex_(columnIndex), rows_(rows) 
{}

template <class Master_matrix>
inline Row_access<Master_matrix>::Row_access(Row_access&& other) noexcept
    : columnIndex_(std::exchange(other.columnIndex_, 0)), rows_(other.rows_) 
{}

template <class Master_matrix>
inline void Row_access<Master_matrix>::insert_entry(ID_index rowIndex, Matrix_entry* entry) 
{
  if (rows_ == nullptr) return;

  if constexpr (!Master_matrix::Option_list::has_removable_rows) {
    if (rows_->size() < rowIndex + 1) rows_->resize(rowIndex + 1);
  }

  // if has_removable_rows should op[] create non existing entry? If not, use try_emplace()
  if constexpr (Master_matrix::Option_list::has_intrusive_rows) {
    rows_->operator[](rowIndex).push_back(*entry);
  } else {
    rows_->operator[](rowIndex).insert(*entry);
  }
}

template <class Master_matrix>
inline void Row_access<Master_matrix>::unlink(Matrix_entry* entry) 
{
  if (rows_ == nullptr) return;

  if constexpr (Master_matrix::Option_list::has_intrusive_rows) {
    entry->Base_hook_matrix_row::unlink();
  } else {
    if constexpr (Master_matrix::Option_list::has_removable_rows) {
      auto it = rows_->find(entry->get_row_index());
      it->second.erase(*entry);
    } else {
      rows_->operator[](entry->get_row_index()).erase(*entry);
    }
  }
}

template <class Master_matrix>
inline void Row_access<Master_matrix>::update_entry(const Matrix_entry& entry) 
{
  if constexpr (!Master_matrix::Option_list::has_intrusive_rows) {
    if (rows_ == nullptr) return;
    auto& row = rows_->at(entry.get_row_index());
    auto it = row.find(entry);
    it = row.erase(it);
    row.insert(it, entry);
  }
}

template <class Master_matrix>
inline typename Row_access<Master_matrix>::Index Row_access<Master_matrix>::get_column_index() const 
{
  return columnIndex_;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_ROW_ACCESS_H
