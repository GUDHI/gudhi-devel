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
 * @file matrix_row_access.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Matrix_row_access class and
 * @ref Gudhi::persistence_matrix::Dummy_matrix_row_access structure.
 */

#ifndef PM_BASE_MATRIX_RA_H
#define PM_BASE_MATRIX_RA_H

#include <utility>  //std::move

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Matrix_row_access, when the the row access is not enabled.
 */
struct Dummy_matrix_row_access
{
  Dummy_matrix_row_access([[maybe_unused]] unsigned int numberOfRows = 0){};

  friend void swap([[maybe_unused]] Dummy_matrix_row_access& d1, [[maybe_unused]] Dummy_matrix_row_access& d2) {}
};

/**
 * @class Matrix_row_access matrix_row_access.h gudhi/Persistence_matrix/matrix_row_access.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the row access for the inheriting matrix.
 * 
 * @tparam Row Either boost::intrusive::list<Matrix_entry,...> if @ref PersistenceMatrixOptions::has_intrusive_rows
 * is true, or std::set<Matrix_entry, RowEntryComp> otherwise.
 * @tparam Row_container Either std::map<Index,Row> if @ref PersistenceMatrixOptions::has_removable_rows is
 *  true, or std::vector<Row> otherwise.
 * @tparam has_removable_rows Value of @ref PersistenceMatrixOptions::has_removable_rows.
 * @tparam ID_index @ref IDIdx index type.
 */
template <typename Row, typename Row_container, bool has_removable_rows, typename ID_index>
class Matrix_row_access 
{
 public:
  /**
   * @brief Default constructor.
   */
  Matrix_row_access() : rows_(new Row_container()){};
  /**
   * @brief Constructor reserving space for the given number of rows.
   *
   * Has only an effect if @ref PersistenceMatrixOptions::has_removable_rows is false.
   * 
   * @param numberOfRows Number of rows to reserve space for.
   */
  Matrix_row_access(unsigned int numberOfRows) : rows_(new Row_container()) {
    if constexpr (!has_removable_rows) {
      rows_->resize(numberOfRows);
    }
  }
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Matrix to copy.
   */
  Matrix_row_access(const Matrix_row_access& toCopy)
      : rows_(new Row_container())  // as the matrix is rebuild, the rows should not be copied.
  {
    if constexpr (!has_removable_rows) {
      rows_->resize(toCopy.rows_->size());
    }
  }
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Matrix_row_access(Matrix_row_access&& other) noexcept : rows_(std::exchange(other.rows_, nullptr)) {}
  /**
   * @brief Destructor.
   */
  ~Matrix_row_access() { delete rows_; }

  /**
   * @brief Returns the row at the given @ref rowindex "row index".
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return: @ref IDIdx for @ref chainmatrix "chain matrices"
   * or updated @ref IDIdx for @ref boundarymatrix "boundary matrices" if swaps occurred.
   * @return Reference to the row.
   */
  Row& get_row(ID_index rowIndex) {
    if constexpr (has_removable_rows) {
      return rows_->at(rowIndex);
    } else {
      return rows_->operator[](rowIndex);
    }
  }
  /**
   * @brief Returns the row at the given @ref rowindex "row index".
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return: @ref IDIdx for @ref chainmatrix "chain matrices"
   * or updated @ref IDIdx for @ref boundarymatrix "boundary matrices" if swaps occurred.
   * @return Const reference to the row.
   */
  const Row& get_row(ID_index rowIndex) const {
    if constexpr (has_removable_rows) {
      return rows_->at(rowIndex);
    } else {
      return rows_->operator[](rowIndex);
    }
  }
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_rows is true. Removes the given row
   * from the row container if the row exists and is empty.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to remove.
   */
  void erase_empty_row(ID_index rowIndex) {
    static_assert(has_removable_rows, "'erase_empty_row' is not implemented for the chosen options.");

    auto it = rows_->find(rowIndex);
    if (it != rows_->end() && it->second.empty()) rows_->erase(it);
  }

  /**
   * @brief Assign operator.
   */
  Matrix_row_access& operator=(const Matrix_row_access& other) {
    if constexpr (has_removable_rows)
      rows_->reserve(other.rows_->size());
    else
      rows_->resize(other.rows_->size());
    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Matrix_row_access& matrix1, Matrix_row_access& matrix2) { std::swap(matrix1.rows_, matrix2.rows_); }

 protected:
  /**
   * @brief Row container. A pointer to facilitate column swaps when two matrices are swapped.
   * Has to be destroyed after matrix_, therefore has to be inherited.
   */
  Row_container* rows_;
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_BASE_MATRIX_RA_H
