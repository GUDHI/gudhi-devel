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
 * @file base_matrix.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Base_matrix class.
 */

#ifndef PM_BASE_MATRIX_H
#define PM_BASE_MATRIX_H

#include <iostream>  //print() only
#include <vector>
#include <utility>  //std::swap, std::move & std::exchange

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Base_matrix base_matrix.h gudhi/Persistence_matrix/base_matrix.h
 * @ingroup persistence_matrix
 *
 * @brief A @ref basematrix "basic matrix" structure allowing to easily manipulate and access entire columns and rows,
 * but not individual cells.
 * 
 * @tparam Master_matrix An instanciation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Base_matrix : public Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >,
                    protected Master_matrix::Matrix_row_access_option 
{
 public:
  using index = typename Master_matrix::index;                        /**< Container index type. */
  using dimension_type = typename Master_matrix::dimension_type;      /**< Dimension value type. */
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element_type = typename Master_matrix::element_type;    /**< Type of a field element. */
  using Column_type = typename Master_matrix::Column_type;            /**< Column type. */
  using container_type = typename Master_matrix::boundary_type;       /**< Type of the column container. */
  using Row_type = typename Master_matrix::Row_type;                  /**< Row type,
                                                                           only necessary with row access option. */
  using Cell_constructor = typename Master_matrix::Cell_constructor;  /**< Factory of @ref Cell classes. */
  using Column_settings = typename Master_matrix::Column_settings;    /**< Structure giving access to the columns to
                                                                           necessary external classes. */

  /**
   * @brief Constructs an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Base_matrix(Column_settings* colSettings);
  /**
   * @brief Constructs a matrix from the given ordered columns. The columns are inserted in the given order.
   * 
   * @tparam Container_type Range type for @ref Matrix::cell_rep_type ranges.
   * Assumed to have a begin(), end() and size() method.
   * @param columns A vector of @ref Matrix::cell_rep_type ranges to construct the columns from.
   * The content of the ranges are assumed to be sorted by increasing ID value.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  template <class Container_type = container_type>
  Base_matrix(const std::vector<Container_type>& columns, 
              Column_settings* colSettings);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Base_matrix(unsigned int numberOfColumns, Column_settings* colSettings);
  /**
   * @brief Copy constructor. If @p colSettings is not a null pointer, its value is kept
   * instead of the one in the copied matrix.
   * 
   * @param matrixToCopy Matrix to copy.
   * @param colSettings Either a pointer to an existing setting structure for the columns or a null pointer.
   * The structure should contain all the necessary external classes specifically necessary for the choosen column type,
   * such as custom allocators. If null pointer, the pointer stored in @p matrixToCopy is used instead.
   */
  Base_matrix(const Base_matrix& matrixToCopy, 
              Column_settings* colSettings = nullptr);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Base_matrix(Base_matrix&& other) noexcept;

  /**
   * @brief Inserts a new ordered column at the end of the matrix by copying the given range of
   * @ref Matrix::cell_rep_type. The content of the range is assumed to be sorted by increasing ID value.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param column Range of @ref Matrix::cell_rep_type from which the column has to be constructed. Assumed to be
   * ordered by increasing ID value. 
   */
  template <class Container_type = container_type>
  void insert_column(const Container_type& column);
  /**
   * @brief Inserts a new ordered column at the given index by copying the given range of @ref Matrix::cell_rep_type.
   * There should not be any other column inserted at that index which was not explicitely removed before.
   * The content of the range is assumed to be sorted by increasing ID value. 
   * Not available when row access is enabled.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param column Range of @ref Matrix::cell_rep_type from which the column has to be constructed. Assumed to be
   * ordered by increasing ID value. 
   * @param columnIndex @ref MatIdx index to which the column has to be inserted.
   */
  template <class Container_type = container_type>
  void insert_column(const Container_type& column, index columnIndex);
  /**
   * @brief Same as @ref insert_column, only for interface purposes. The given dimension is ignored and not stored.
   * 
   * @tparam Boundary_type Range of @ref Matrix::cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param boundary Range of @ref Matrix::cell_rep_type from which the column has to be constructed. Assumed to be
   * ordered by increasing ID value. 
   * @param dim Ignored.
   */
  template <class Boundary_type>
  void insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
  /**
   * @brief Returns the column at the given @ref MatIdx index.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   *
   * Note that before returning the column, all column cells can eventually be reordered, if lazy swaps occurred.
   * It is therefore recommended to avoid calling @ref get_column between column or row swaps, otherwise the benefits
   * of the the lazyness is lost.
   * 
   * @param columnIndex @ref MatIdx index of the column to return.
   * @return Reference to the column.
   */
  Column_type& get_column(index columnIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access is true.
   * Returns the row at the given @ref rowindex "row index" of the matrix.
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   *
   * Note that before returning the row, all column cells can eventually be reordered, if lazy swaps occurred.
   * It is therefore recommended to avoid calling @ref get_row between column or row swaps, otherwise the benefits
   * of the the lazyness is lost.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return.
   * @return Reference to the row.
   */
  Row_type& get_row(index rowIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_map_column_container is true. Otherwise, see
   * @ref remove_last. Erases the given column from the matrix. If the given column index corresponded to the last
   * used column index, the "new last column index" will be `columnIndex - 1`, even if a column was never explicitely
   * inserted at this index (possible when
   * @ref insert_column(const Container_type& column, index columnIndex) "insert_column(column, columnIndex)" was used).
   * If the column didn't existed, it will simply be considered as an empty column.
   *
   * If @ref PersistenceMatrixOptions::has_row_access is also true, the deleted column cells are also automatically
   * removed from their  respective rows.
   * 
   * @param columnIndex @ref MatIdx index of the column to remove.
   */
  void remove_column(index columnIndex);
  /**
   * @brief Removes the last column from the matrix. The last column is at index \f$ max(ins1, ins2, rem - 1) \f$,
   * where:
   * - \f$ ins1 \f$ is the index of the last column inserted by
   * @ref insert_column(const Container_type& column) "insert_column(column)",
   * - \f$ ins2 \f$ is largest index used with
   * @ref insert_column(const Container_type& column, index columnIndex) "insert_column(column, columnIndex)",
   * - \f$ rem \f$ is the last index just before the last call of @ref remove_column or @ref remove_last.
   *
   * If \f$ max(ins1, ins2, rem - 1) = rem - 1 \f$ but no column was explicitely inserted at that index (possible
   * by the use of
   * @ref insert_column(const Container_type& column, index columnIndex) "insert_column(column, columnIndex)"),
   * the column is assumed to be an empty column.
   *
   * See also @ref remove_column.
   */
  void remove_last();
  /**
   * @brief If @ref PersistenceMatrixOptions::has_row_access and @ref PersistenceMatrixOptions::has_removable_rows
   * are true: assumes that the row is empty and removes it. If @ref PersistenceMatrixOptions::has_map_column_container
   * and @ref PersistenceMatrixOptions::has_column_and_row_swaps are true: cleans up maps used for the lazy row swaps.
   * Otherwise, does nothing.
   *
   * @warning The removed rows are always assumed to be empty. If it is not the case, the deleted row cells are not
   * removed from their columns. And in the case of intrusive rows, this will generate a segmentation fault when 
   * the column cells are destroyed later. The row access is just meant as a "read only" access to the rows and the
   * @ref erase_empty_row method just as a way to specify that a row is empty and can therefore be removed from
   * dictionnaries. This allows to avoid testing the emptiness of a row at each column cell removal, what can be
   * quite frequent. 
   * 
   * @param rowIndex @ref rowindex "Row index" of the empty row.
   */
  void erase_empty_row(index rowIndex);

  /**
   * @brief Returns the current number of columns in the matrix.
   * 
   * @return The number of columns.
   */
  index get_number_of_columns() const;

  /**
   * @brief Adds column represented by @p sourceColumn onto the column at @p targetColumnIndex in the matrix.
   * 
   * @tparam Cell_range_or_column_index Either a range of @ref Cell with a begin() and end() method,
   * or any integer type.
   * @param sourceColumn Either a cell range or the @ref MatIdx index of the column to add.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Cell_range_or_column_index>
  void add_to(const Cell_range_or_column_index& sourceColumn, index targetColumnIndex);
  /**
   * @brief Multiplies the target column with the coefficiant and then adds the source column to it.
   * That is: `targetColumn = (targetColumn * coefficient) + sourceColumn`.
   * 
   * @tparam Cell_range_or_column_index Either a range of @ref Cell with a begin() and end() method,
   * or any integer type.
   * @param sourceColumn Either a cell range or the @ref MatIdx index of the column to add.
   * @param coefficient Value to multiply.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Cell_range_or_column_index>
  void multiply_target_and_add_to(const Cell_range_or_column_index& sourceColumn, 
                                  const Field_element_type& coefficient,
                                  index targetColumnIndex);
  /**
   * @brief Multiplies the source column with the coefficiant before adding it to the target column.
   * That is: `targetColumn += (coefficient * sourceColumn)`. The source column will **not** be modified.
   * 
   * @tparam Cell_range_or_column_index Either a range of @ref Cell with a begin() and end() method,
   * or any integer type.
   * @param coefficient Value to multiply.
   * @param sourceColumn Either a cell range or the @ref MatIdx index of the column to add.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Cell_range_or_column_index>
  void multiply_source_and_add_to(const Field_element_type& coefficient, 
                                  const Cell_range_or_column_index& sourceColumn,
                                  index targetColumnIndex);

  /**
   * @brief Zeroes the cell at the given coordinates.
   * 
   * @param columnIndex @ref MatIdx index of the column of the cell.
   * @param rowIndex @ref rowindex "Row index" of the row of the cell.
   */
  void zero_cell(index columnIndex, index rowIndex);
  /**
   * @brief Zeroes the column at the given index.
   * 
   * @param columnIndex @ref MatIdx index of the column to zero.
   */
  void zero_column(index columnIndex);
  /**
   * @brief Indicates if the cell at given coordinates has value zero.
   * 
   * @param columnIndex @ref MatIdx index of the column of the cell.
   * @param rowIndex @ref rowindex "Row index" of the row of the cell.
   * @return true If the cell has value zero.
   * @return false Otherwise.
   */
  bool is_zero_cell(index columnIndex, index rowIndex) const;
  /**
   * @brief Indicates if the column at given index has value zero.
   * 
   * @param columnIndex @ref MatIdx index of the column.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(index columnIndex);

  /**
   * @brief Resets the matrix to an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  void reset(Column_settings* colSettings) {
    matrix_.clear();
    nextInsertIndex_ = 0;
    colSettings_ = colSettings;
  }

  /**
   * @brief Assign operator.
   */
  Base_matrix& operator=(const Base_matrix& other);
  /**
   * @brief Swap operator.
   */
  friend void swap(Base_matrix& matrix1, Base_matrix& matrix2) {
    swap(static_cast<typename Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >&>(matrix1),
         static_cast<typename Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >&>(matrix2));
    matrix1.matrix_.swap(matrix2.matrix_);
    std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
    std::swap(matrix1.colSettings_, matrix2.colSettings_);

    if constexpr (Master_matrix::Option_list::has_row_access) {
      swap(static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix1),
           static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix2));
    }
  }

  void print();  // for debug

 private:
  using swap_opt = typename Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >;
  using ra_opt = typename Master_matrix::Matrix_row_access_option;
  using matrix_type = typename Master_matrix::column_container_type;
  using cell_rep_type =
      typename std::conditional<Master_matrix::Option_list::is_z2, 
                                index, 
                                std::pair<index, Field_element_type>
                               >::type;

  friend swap_opt;  // direct access to matrix_ to avoid row reorder.

  matrix_type matrix_;            /**< Column container. */
  index nextInsertIndex_;         /**< Next unused column index. */
  Column_settings* colSettings_;  /**< Cell factory. */

  template <class Container_type = container_type>
  void _insert(const Container_type& column, index columnIndex, dimension_type dim);
  void _orderRowsIfNecessary();
  const Column_type& _get_column(index columnIndex) const;
  Column_type& _get_column(index columnIndex);
  index _get_real_row_index(index rowIndex) const;
  template <class Container_type>
  void _container_insert(const Container_type& column, index pos, dimension_type dim);
  void _container_insert(const Column_type& column, [[maybe_unused]] index pos = 0);
};

template <class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(Column_settings* colSettings)
    : swap_opt(), ra_opt(), nextInsertIndex_(0), colSettings_(colSettings)
{}

template <class Master_matrix>
template <class Container_type>
inline Base_matrix<Master_matrix>::Base_matrix(const std::vector<Container_type>& columns, 
                                               Column_settings* colSettings)
    : swap_opt(columns.size()),
      // not ideal if max row index is much smaller than max column index, does that happen often?
      ra_opt(columns.size()),
      matrix_(!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access
                  ? 0
                  : columns.size()),
      nextInsertIndex_(columns.size()),
      colSettings_(colSettings)
{
  if constexpr (!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access)
    matrix_.reserve(columns.size());

  for (index i = 0; i < columns.size(); i++) {
    _container_insert(columns[i], i, columns[i].size() == 0 ? 0 : columns[i].size() - 1);
  }
}

template <class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(unsigned int numberOfColumns, 
                                               Column_settings* colSettings)
    : swap_opt(numberOfColumns),
      ra_opt(numberOfColumns),
      matrix_(!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access
                  ? 0
                  : numberOfColumns),
      nextInsertIndex_(0),
      colSettings_(colSettings)
{
  if constexpr (!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access)
    matrix_.reserve(numberOfColumns);
}

template <class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(const Base_matrix& matrixToCopy, 
                                               Column_settings* colSettings)
    : swap_opt(static_cast<const swap_opt&>(matrixToCopy)),
      ra_opt(static_cast<const ra_opt&>(matrixToCopy)),
      nextInsertIndex_(matrixToCopy.nextInsertIndex_),
      colSettings_(colSettings == nullptr ? matrixToCopy.colSettings_ : colSettings)
{
  matrix_.reserve(matrixToCopy.matrix_.size());
  for (const auto& cont : matrixToCopy.matrix_){
    if constexpr (Master_matrix::Option_list::has_map_column_container){
      _container_insert(cont.second, cont.first);
    } else {
      _container_insert(cont);
    }
  }
}

template <class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(Base_matrix&& other) noexcept
    : swap_opt(std::move(static_cast<swap_opt&>(other))),
      ra_opt(std::move(static_cast<ra_opt&>(other))),
      matrix_(std::move(other.matrix_)),
      nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0)),
      colSettings_(std::exchange(other.colSettings_, nullptr)) 
{}

template <class Master_matrix>
template <class Container_type>
inline void Base_matrix<Master_matrix>::insert_column(const Container_type& column) 
{
  //TODO: dim not actually stored right now, so either get rid of it or store it again
  _insert(column, nextInsertIndex_, column.size() == 0 ? 0 : column.size() - 1);
  ++nextInsertIndex_;
}

template <class Master_matrix>
template <class Container_type>
inline void Base_matrix<Master_matrix>::insert_column(const Container_type& column, index columnIndex) 
{
  static_assert(!Master_matrix::Option_list::has_row_access,
                "Columns have to be inserted at the end of the matrix when row access is enabled.");

  if (columnIndex >= nextInsertIndex_) nextInsertIndex_ = columnIndex + 1;
  //TODO: dim not actually stored right now, so either get rid of it or store it again
  _insert(column, columnIndex, column.size() == 0 ? 0 : column.size() - 1);
}

template <class Master_matrix>
template <class Boundary_type>
inline void Base_matrix<Master_matrix>::insert_boundary(const Boundary_type& boundary, dimension_type dim) 
{
  if (dim == -1) dim = boundary.size() == 0 ? 0 : boundary.size() - 1;
  //TODO: dim not actually stored right now, so either get rid of it or store it again
  _insert(boundary, nextInsertIndex_++, dim);
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::Column_type& Base_matrix<Master_matrix>::get_column(index columnIndex) 
{
  _orderRowsIfNecessary();
  return _get_column(columnIndex);
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::Row_type& Base_matrix<Master_matrix>::get_row(index rowIndex) 
{
  static_assert(Master_matrix::Option_list::has_row_access, "Row access has to be enabled for this method.");

  _orderRowsIfNecessary();
  return ra_opt::get_row(rowIndex);
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::remove_column(index columnIndex) 
{
  static_assert(Master_matrix::Option_list::has_map_column_container,
                "'remove_column' is not implemented for the chosen options.");

  // assumes that eventual "holes" left at unused indices are considered as empty columns.
  if (columnIndex == nextInsertIndex_ - 1) --nextInsertIndex_;

  matrix_.erase(columnIndex);
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::remove_last() 
{
  if (nextInsertIndex_ == 0) return;  //empty matrix
  --nextInsertIndex_;  // assumes that eventual "holes" left at unused indices are considered as empty columns.

  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    matrix_.erase(nextInsertIndex_);
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      GUDHI_CHECK(nextInsertIndex_ == matrix_.size() - 1,
                  std::logic_error("Base_matrix::remove_last - Indexation problem."));
      matrix_.pop_back();
    } else {
      matrix_[nextInsertIndex_].clear();
    }
  }
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::erase_empty_row(index rowIndex) 
{
  if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows) {
    ra_opt::erase_empty_row(_get_real_row_index(rowIndex));
  }
  if constexpr ((Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) &&
                Master_matrix::Option_list::has_map_column_container) {
    auto it = swap_opt::indexToRow_.find(rowIndex);
    swap_opt::rowToIndex_.erase(it->second);
    swap_opt::indexToRow_.erase(it);
  }
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::index Base_matrix<Master_matrix>::get_number_of_columns() const 
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.size();
  } else {
    return nextInsertIndex_;  // matrix could have been resized much bigger while insert
  }
}

template <class Master_matrix>
template <class Cell_range_or_column_index>
inline void Base_matrix<Master_matrix>::add_to(const Cell_range_or_column_index& sourceColumn,
                                               index targetColumnIndex) 
{
  if constexpr (std::is_integral_v<Cell_range_or_column_index>) {
    _get_column(targetColumnIndex) += _get_column(sourceColumn);
  } else {
    _get_column(targetColumnIndex) += sourceColumn;
  }
}

template <class Master_matrix>
template <class Cell_range_or_column_index>
inline void Base_matrix<Master_matrix>::multiply_target_and_add_to(const Cell_range_or_column_index& sourceColumn,
                                                                   const Field_element_type& coefficient,
                                                                   index targetColumnIndex) 
{
  if constexpr (std::is_integral_v<Cell_range_or_column_index>) {
    _get_column(targetColumnIndex).multiply_target_and_add(coefficient, _get_column(sourceColumn));
  } else {
    _get_column(targetColumnIndex).multiply_target_and_add(coefficient, sourceColumn);
  }
}

template <class Master_matrix>
template <class Cell_range_or_column_index>
inline void Base_matrix<Master_matrix>::multiply_source_and_add_to(const Field_element_type& coefficient,
                                                                   const Cell_range_or_column_index& sourceColumn,
                                                                   index targetColumnIndex) 
{
  if constexpr (std::is_integral_v<Cell_range_or_column_index>) {
    _get_column(targetColumnIndex).multiply_source_and_add(_get_column(sourceColumn), coefficient);
  } else {
    _get_column(targetColumnIndex).multiply_source_and_add(sourceColumn, coefficient);
  }
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_cell(index columnIndex, index rowIndex) 
{
  _get_column(columnIndex).clear(_get_real_row_index(rowIndex));
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_column(index columnIndex) {
  _get_column(columnIndex).clear();
}

template <class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const 
{
  return !(_get_column(columnIndex).is_non_zero(_get_real_row_index(rowIndex)));
}

template <class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_column(index columnIndex) 
{
  return _get_column(columnIndex).is_empty();
}

template <class Master_matrix>
inline Base_matrix<Master_matrix>& Base_matrix<Master_matrix>::operator=(const Base_matrix& other) 
{
  swap_opt::operator=(other);
  ra_opt::operator=(other);
  matrix_.clear();
  nextInsertIndex_ = other.nextInsertIndex_;
  colSettings_ = other.colSettings_;

  matrix_.reserve(other.matrix_.size());
  for (const auto& cont : other.matrix_){
    if constexpr (Master_matrix::Option_list::has_map_column_container){
      _container_insert(cont.second, cont.first);
    } else {
      _container_insert(cont);
    }
  }

  return *this;
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::print() 
{
  _orderRowsIfNecessary();
  std::cout << "Base_matrix:\n";
  for (index i = 0; i < nextInsertIndex_; ++i) {
    const Column_type& col = matrix_[i];
    for (const auto& e : col.get_content(nextInsertIndex_)) {
      if (e == 0u)
        std::cout << "- ";
      else
        std::cout << e << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  if constexpr (Master_matrix::Option_list::has_row_access) {
    std::cout << "Row Matrix:\n";
    for (index i = 0; i < nextInsertIndex_; ++i) {
      const auto& row = ra_opt::rows_[i];
      for (const auto& cell : row) {
        std::cout << cell.get_column_index() << " ";
      }
      std::cout << "(" << i << ")\n";
    }
    std::cout << "\n";
  }
}

template <class Master_matrix>
template <class Container_type>
inline void Base_matrix<Master_matrix>::_insert(const Container_type& column, index columnIndex, dimension_type dim) 
{
  _orderRowsIfNecessary();

  //resize of containers when necessary:
  index pivot = 0;
  if (column.begin() != column.end()) {
    //first, compute pivot of `column`
    if constexpr (Master_matrix::Option_list::is_z2) {
      pivot = *std::prev(column.end());
    } else {
      pivot = std::prev(column.end())->first;
    }
    //row container
    if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_removable_rows)
      if (ra_opt::rows_->size() <= pivot) ra_opt::rows_->resize(pivot + 1);
  }

  //row swap map containers
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
      for (auto id : column) {
        index idx;
        if constexpr (Master_matrix::Option_list::is_z2) {
          idx = id;
        } else {
          idx = id.first;
        }
        swap_opt::indexToRow_[idx] = idx;
        swap_opt::rowToIndex_[idx] = idx;
      }
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
      index size = swap_opt::indexToRow_.size();
      if (size <= pivot) {
        for (index i = size; i <= pivot; i++) {
          swap_opt::indexToRow_.push_back(i);
          swap_opt::rowToIndex_.push_back(i);
        }
      }
    }
    //column container
    if constexpr (!Master_matrix::Option_list::has_row_access) {
      if (matrix_.size() <= columnIndex) {
        matrix_.resize(columnIndex + 1);
      }
    }
  }

  _container_insert(column, columnIndex, dim);
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::_orderRowsIfNecessary() 
{
  if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
    if (swap_opt::rowSwapped_) swap_opt::_orderRows();
  }
}

template <class Master_matrix>
inline const typename Base_matrix<Master_matrix>::Column_type& Base_matrix<Master_matrix>::_get_column(
    index columnIndex) const
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.at(columnIndex);
  } else {
    return matrix_[columnIndex];
  }
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::Column_type& Base_matrix<Master_matrix>::_get_column(index columnIndex)
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.at(columnIndex);
  } else {
    return matrix_[columnIndex];
  }
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::index Base_matrix<Master_matrix>::_get_real_row_index(index rowIndex) const
{
  if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      return swap_opt::indexToRow_.at(rowIndex);
    } else {
      return swap_opt::indexToRow_[rowIndex];
    }
  } else {
    return rowIndex;
  }
}

template <class Master_matrix>
template <class Container_type>
inline void Base_matrix<Master_matrix>::_container_insert(const Container_type& column, index pos, dimension_type dim){
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.try_emplace(pos, Column_type(pos, column, dim, ra_opt::rows_, colSettings_));
    } else {
      matrix_.try_emplace(pos, Column_type(column, dim, colSettings_));
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.emplace_back(pos, column, dim, ra_opt::rows_, colSettings_);
    } else {
      matrix_[pos] = Column_type(column, dim, colSettings_);
    }
  }
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::_container_insert(const Column_type& column, [[maybe_unused]] index pos){
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.try_emplace(pos, Column_type(column, column.get_column_index(), ra_opt::rows_, colSettings_));
    } else {
      matrix_.try_emplace(pos, Column_type(column, colSettings_));
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.emplace_back(column, column.get_column_index(), ra_opt::rows_, colSettings_);
    } else {
      matrix_.emplace_back(column, colSettings_);
    }
  }
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_BASE_MATRIX_H
