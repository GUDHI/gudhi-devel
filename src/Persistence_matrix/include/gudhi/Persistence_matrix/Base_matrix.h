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
 * @file Base_matrix.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Base_matrix class.
 */

#ifndef PM_BASE_MATRIX_H
#define PM_BASE_MATRIX_H

#include <iostream>  //print() only
#include <vector>
#include <utility>  //std::swap, std::move & std::exchange

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Base_matrix Base_matrix.h gudhi/Persistence_matrix/Base_matrix.h
 * @ingroup persistence_matrix
 *
 * @brief A @ref basematrix "basic matrix" structure allowing to easily manipulate and access entire columns and rows,
 * but not individual entries.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Base_matrix : public Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >,
                    protected Master_matrix::Matrix_row_access_option
{
 private:
  using Swap_opt = typename Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >;
  using RA_opt = typename Master_matrix::Matrix_row_access_option;

 public:
  using Index = typename Master_matrix::Index;         /**< Container index type. */
  using Dimension = typename Master_matrix::Dimension; /**< Dimension value type. */
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element = typename Master_matrix::Element;               /**< Type of a field element. */
  using Column = typename Master_matrix::Column;                       /**< Column type. */
  using Boundary = typename Master_matrix::Boundary;                   /**< Type of the column container. */
  using Row = typename Master_matrix::Row;                             /**< Row type,
                                                                            only necessary with row access option. */
  using Entry_constructor = typename Master_matrix::Entry_constructor; /**< Factory of @ref Entry classes. */
  using Column_settings = typename Master_matrix::Column_settings;     /**< Structure giving access to the columns to
                                                                            necessary external classes. */

  /**
   * @brief Constructs an empty matrix.
   *
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  Base_matrix(Column_settings* colSettings);
  /**
   * @brief Constructs a matrix from the given ordered columns. The columns are inserted in the given order.
   *
   * @tparam Container Range type for @ref Matrix::Entry_representative ranges.
   * Assumed to have a begin(), end() and size() method.
   * @param columns A vector of @ref Matrix::Entry_representative ranges to construct the columns from.
   * The content of the ranges are assumed to be sorted by increasing ID value.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  template <class Container = Boundary>
  Base_matrix(const std::vector<Container>& columns, Column_settings* colSettings);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   *
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  Base_matrix(unsigned int numberOfColumns, Column_settings* colSettings);
  /**
   * @brief Copy constructor. If @p colSettings is not a null pointer, its value is kept
   * instead of the one in the copied matrix.
   *
   * @param matrixToCopy Matrix to copy.
   * @param colSettings Either a pointer to an existing setting structure for the columns or a null pointer.
   * The structure should contain all the necessary external classes specifically necessary for the chosen column type,
   * such as custom allocators. If null pointer, the pointer stored in @p matrixToCopy is used instead.
   */
  Base_matrix(const Base_matrix& matrixToCopy, Column_settings* colSettings = nullptr);
  /**
   * @brief Move constructor.
   *
   * @param other Matrix to move.
   */
  Base_matrix(Base_matrix&& other) noexcept;

  ~Base_matrix() = default;

  /**
   * @brief Inserts a new ordered column at the end of the matrix by copying the given range of
   * @ref Matrix::Entry_representative. The content of the range is assumed to be sorted by increasing ID value.
   *
   * @tparam Container Range of @ref Matrix::Entry_representative. Assumed to have a begin(), end() and size() method.
   * @param column Range of @ref Matrix::Entry_representative from which the column has to be constructed. Assumed to be
   * ordered by increasing ID value.
   */
  template <class Container = Boundary, class = std::enable_if_t<!std::is_arithmetic_v<Container> > >
  void insert_column(const Container& column);
  /**
   * @brief Inserts a new ordered column at the given index by copying the given range of
   * @ref Matrix::Entry_representative.
   * There should not be any other column inserted at that index which was not explicitly removed before.
   * The content of the range is assumed to be sorted by increasing ID value.
   * Not available when row access is enabled.
   *
   * @tparam Container Range of @ref Matrix::Entry_representative. Assumed to have a begin(), end() and size() method.
   * @param column Range of @ref Matrix::Entry_representative from which the column has to be constructed. Assumed to be
   * ordered by increasing ID value.
   * @param columnIndex @ref MatIdx index to which the column has to be inserted.
   */
  template <class Container = Boundary, class = std::enable_if_t<!std::is_arithmetic_v<Container> > >
  void insert_column(const Container& column, Index columnIndex);
  /**
   * @brief Inserts a new column at the end of the matrix. The column will consist of the given index only.
   * 
   * @param idx Entry ID.
   * @param e Entry coefficient. Ignored if the coefficient field is Z2. Default value: 0.
   */
  void insert_column(Index idx, Field_element e = 0);
  /**
   * @brief Same as @ref insert_column, only for interface purposes. The given dimension is ignored and not stored.
   *
   * @tparam Boundary_range Range of @ref Matrix::Entry_representative. Assumed to have a begin(), end() and size()
   * method.
   * @param boundary Range of @ref Matrix::Entry_representative from which the column has to be constructed. Assumed to
   * be ordered by increasing ID value.
   * @param dim Ignored.
   */
  template <class Boundary_range>
  void insert_boundary(const Boundary_range& boundary,
                       Dimension dim = Master_matrix::template get_null_value<Dimension>());
  /**
   * @brief Returns the column at the given @ref MatIdx index.
   * The type of the column depends on the chosen options, see @ref PersistenceMatrixOptions::column_type.
   *
   * Note that before returning the column, all column entries can eventually be reordered, if lazy swaps occurred.
   * It is therefore recommended to avoid calling @ref get_column between column or row swaps, otherwise the benefits
   * of the the laziness is lost.
   *
   * @param columnIndex @ref MatIdx index of the column to return.
   * @return Reference to the column.
   */
  Column& get_column(Index columnIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access is true.
   * Returns the row at the given @ref rowindex "row index" of the matrix.
   * The type of the row depends on the chosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   *
   * Note that before returning the row, all column entries can eventually be reordered, if lazy swaps occurred.
   * It is therefore recommended to avoid calling @ref get_row between column or row swaps, otherwise the benefits
   * of the the laziness is lost.
   *
   * @param rowIndex @ref rowindex "Row index" of the row to return.
   * @return Reference to the row.
   */
  Row& get_row(Index rowIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_map_column_container is true. Otherwise, see
   * @ref remove_last. Erases the given column from the matrix. If the given column index corresponded to the last
   * used column index, the "new last column index" will be `columnIndex - 1`, even if a column was never explicitly
   * inserted at this index (possible when
   * @ref insert_column(const Container& column, Index columnIndex) "insert_column(column, columnIndex)" was used).
   * If the column didn't existed, it will simply be considered as an empty column.
   *
   * If @ref PersistenceMatrixOptions::has_row_access is also true, the deleted column entries are also automatically
   * removed from their  respective rows.
   *
   * @param columnIndex @ref MatIdx index of the column to remove.
   */
  void remove_column(Index columnIndex);
  /**
   * @brief Removes the last column from the matrix. The last column is at index \f$ max(ins1, ins2, rem - 1) \f$,
   * where:
   * - \f$ ins1 \f$ is the index of the last column inserted by
   * @ref insert_column(const Container& column) "insert_column(column)",
   * - \f$ ins2 \f$ is largest index used with
   * @ref insert_column(const Container& column, Index columnIndex) "insert_column(column, columnIndex)",
   * - \f$ rem \f$ is the last index just before the last call of @ref remove_column or @ref remove_last.
   *
   * If \f$ max(ins1, ins2, rem - 1) = rem - 1 \f$ but no column was explicitly inserted at that index (possible
   * by the use of
   * @ref insert_column(const Container& column, Index columnIndex) "insert_column(column, columnIndex)"),
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
   * @warning The removed rows are always assumed to be empty. If it is not the case, the deleted row entries are not
   * removed from their columns. And in the case of intrusive rows, this will generate a segmentation fault when
   * the column entries are destroyed later. The row access is just meant as a "read only" access to the rows and the
   * @ref erase_empty_row method just as a way to specify that a row is empty and can therefore be removed from
   * dictionaries. This allows to avoid testing the emptiness of a row at each column entry removal, what can be
   * quite frequent.
   *
   * @param rowIndex @ref rowindex "Row index" of the empty row.
   */
  void erase_empty_row(Index rowIndex);

  /**
   * @brief Returns the current number of columns in the matrix.
   *
   * @return The number of columns.
   */
  Index get_number_of_columns() const;

  /**
   * @brief Adds column represented by @p sourceColumn onto the column at @p targetColumnIndex in the matrix.
   *
   * @tparam Entry_range_or_column_index Either a range of @ref Entry with a begin() and end() method,
   * or any integer type.
   * @param sourceColumn Either an entry range or the @ref MatIdx index of the column to add.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Entry_range_or_column_index>
  void add_to(const Entry_range_or_column_index& sourceColumn, Index targetColumnIndex);
  /**
   * @brief Multiplies the target column with the coefficient and then adds the source column to it.
   * That is: `targetColumn = (targetColumn * coefficient) + sourceColumn`.
   *
   * @tparam Entry_range_or_column_index Either a range of @ref Entry with a begin() and end() method,
   * or any integer type.
   * @param sourceColumn Either an entry range or the @ref MatIdx index of the column to add.
   * @param coefficient Value to multiply.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Entry_range_or_column_index>
  void multiply_target_and_add_to(const Entry_range_or_column_index& sourceColumn,
                                  const Field_element& coefficient,
                                  Index targetColumnIndex);
  /**
   * @brief Multiplies the source column with the coefficient before adding it to the target column.
   * That is: `targetColumn += (coefficient * sourceColumn)`. The source column will **not** be modified.
   *
   * @tparam Entry_range_or_column_index Either a range of @ref Entry with a begin() and end() method,
   * or any integer type.
   * @param coefficient Value to multiply.
   * @param sourceColumn Either an entry range or the @ref MatIdx index of the column to add.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Entry_range_or_column_index>
  void multiply_source_and_add_to(const Field_element& coefficient,
                                  const Entry_range_or_column_index& sourceColumn,
                                  Index targetColumnIndex);

  /**
   * @brief Zeroes the entry at the given coordinates.
   *
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   */
  void zero_entry(Index columnIndex, Index rowIndex);
  /**
   * @brief Zeroes the column at the given index.
   *
   * @param columnIndex @ref MatIdx index of the column to zero.
   */
  void zero_column(Index columnIndex);
  /**
   * @brief Indicates if the entry at given coordinates has value zero.
   *
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @return true If the entry has value zero.
   * @return false Otherwise.
   */
  bool is_zero_entry(Index columnIndex, Index rowIndex) const;
  /**
   * @brief Indicates if the column at given index has value zero.
   *
   * @param columnIndex @ref MatIdx index of the column.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(Index columnIndex);

  /**
   * @brief Resets the matrix to an empty matrix.
   *
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  void reset(Column_settings* colSettings)
  {
    if constexpr (Master_matrix::Option_list::has_vine_update || Master_matrix::Option_list::has_column_and_row_swaps)
      Swap_opt::_reset();
    matrix_.clear();
    nextInsertIndex_ = 0;
    colSettings_ = colSettings;
  }

  /**
   * @brief Assign operator.
   */
  Base_matrix& operator=(const Base_matrix& other);
  /**
   * @brief Move assign operator.
   */
  Base_matrix& operator=(Base_matrix&& other) noexcept;

  /**
   * @brief Swap operator.
   */
  friend void swap(Base_matrix& matrix1, Base_matrix& matrix2) noexcept
  {
    swap(static_cast<Swap_opt&>(matrix1), static_cast<Swap_opt&>(matrix2));
    matrix1.matrix_.swap(matrix2.matrix_);
    std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
    std::swap(matrix1.colSettings_, matrix2.colSettings_);

    if constexpr (Master_matrix::Option_list::has_row_access) {
      swap(static_cast<RA_opt&>(matrix1), static_cast<RA_opt&>(matrix2));
    }
  }

  void print();  // for debug

 private:
  using Column_container = typename Master_matrix::Column_container;
  using Entry_representative =
      typename std::conditional<Master_matrix::Option_list::is_z2, Index, std::pair<Index, Field_element> >::type;

  friend Swap_opt;  // direct access to matrix_ to avoid row reorder.

  Column_container matrix_;      /**< Column container. */
  Index nextInsertIndex_;        /**< Next unused column index. */
  Column_settings* colSettings_; /**< Entry factory. */

  template <class Container = Boundary, class = std::enable_if_t<!std::is_arithmetic_v<Container> > >
  void _insert(const Container& column, Index columnIndex, Dimension dim);
  void _insert(Index idx, Field_element e, Index columnIndex, Dimension dim);
  void _orderRowsIfNecessary();
  const Column& _get_column(Index columnIndex) const;
  Column& _get_column(Index columnIndex);
  Index _get_real_row_index(Index rowIndex) const;
  template <class Container>
  void _container_insert(const Container& column, Index pos, Dimension dim);
  void _container_insert(Index idx, Field_element e, Index pos, Dimension dim);
  template <class ColumnIterator>
  void _container_insert(const ColumnIterator& rep);
};

template <class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(Column_settings* colSettings)
    : Swap_opt(), RA_opt(), nextInsertIndex_(0), colSettings_(colSettings)
{}

template <class Master_matrix>
template <class Container>
inline Base_matrix<Master_matrix>::Base_matrix(const std::vector<Container>& columns, Column_settings* colSettings)
    : Swap_opt(columns.size()),
      // not ideal if max row index is much smaller than max column index, does that happen often?
      RA_opt(columns.size()),
      matrix_(!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access
                  ? 0
                  : columns.size()),
      nextInsertIndex_(columns.size()),
      colSettings_(colSettings)
{
  if constexpr (!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access)
    matrix_.reserve(columns.size());

  for (Index i = 0; i < columns.size(); i++) {
    _container_insert(columns[i], i, columns[i].size() == 0 ? 0 : columns[i].size() - 1);
  }
}

template <class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(unsigned int numberOfColumns, Column_settings* colSettings)
    : Swap_opt(numberOfColumns),
      RA_opt(numberOfColumns),
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
inline Base_matrix<Master_matrix>::Base_matrix(const Base_matrix& matrixToCopy, Column_settings* colSettings)
    : Swap_opt(static_cast<const Swap_opt&>(matrixToCopy)),
      RA_opt(static_cast<const RA_opt&>(matrixToCopy)),
      nextInsertIndex_(matrixToCopy.nextInsertIndex_),
      colSettings_(colSettings == nullptr ? matrixToCopy.colSettings_ : colSettings)
{
  matrix_.reserve(matrixToCopy.matrix_.size());
  for (const auto& cont : matrixToCopy.matrix_) {
    _container_insert(cont);
  }
}

template <class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(Base_matrix&& other) noexcept
    : Swap_opt(std::move(static_cast<Swap_opt&>(other))),
      RA_opt(std::move(static_cast<RA_opt&>(other))),
      matrix_(std::move(other.matrix_)),
      nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0)),
      colSettings_(std::exchange(other.colSettings_, nullptr))
{
}

template <class Master_matrix>
template <class Container, class>
inline void Base_matrix<Master_matrix>::insert_column(const Container& column)
{
  // TODO: dim not actually stored right now, so either get rid of it or store it again
  _insert(column, nextInsertIndex_, column.size() == 0 ? 0 : column.size() - 1);
  ++nextInsertIndex_;
}

template <class Master_matrix>
template <class Container, class>
inline void Base_matrix<Master_matrix>::insert_column(const Container& column, Index columnIndex)
{
  static_assert(!Master_matrix::Option_list::has_row_access,
                "Columns have to be inserted at the end of the matrix when row access is enabled.");

  if (columnIndex >= nextInsertIndex_) nextInsertIndex_ = columnIndex + 1;
  // TODO: dim not actually stored right now, so either get rid of it or store it again
  _insert(column, columnIndex, column.size() == 0 ? 0 : column.size() - 1);
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::insert_column(Index idx, Field_element e)
{
  // TODO: dim not actually stored right now, so either get rid of it or store it again
  _insert(idx, e, nextInsertIndex_, 0);
  ++nextInsertIndex_;
}

template <class Master_matrix>
template <class Boundary_range>
inline void Base_matrix<Master_matrix>::insert_boundary(const Boundary_range& boundary, Dimension dim)
{
  if (dim == Master_matrix::template get_null_value<Dimension>()) dim = boundary.size() == 0 ? 0 : boundary.size() - 1;
  // TODO: dim not actually stored right now, so either get rid of it or store it again
  _insert(boundary, nextInsertIndex_++, dim);
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::Column& Base_matrix<Master_matrix>::get_column(Index columnIndex)
{
  _orderRowsIfNecessary();
  return _get_column(columnIndex);
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::Row& Base_matrix<Master_matrix>::get_row(Index rowIndex)
{
  static_assert(Master_matrix::Option_list::has_row_access, "Row access has to be enabled for this method.");

  _orderRowsIfNecessary();
  return RA_opt::get_row(rowIndex);
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::remove_column(Index columnIndex)
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
  if (nextInsertIndex_ == 0) return;  // empty matrix
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
inline void Base_matrix<Master_matrix>::erase_empty_row(Index rowIndex)
{
  if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows) {
    RA_opt::erase_empty_row(_get_real_row_index(rowIndex));
  }
  if constexpr ((Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) &&
                Master_matrix::Option_list::has_map_column_container) {
    Swap_opt::_erase_row(rowIndex);
  }
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::Index Base_matrix<Master_matrix>::get_number_of_columns() const
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.size();
  } else {
    return nextInsertIndex_;  // matrix could have been resized much bigger while insert
  }
}

template <class Master_matrix>
template <class Entry_range_or_column_index>
inline void Base_matrix<Master_matrix>::add_to(const Entry_range_or_column_index& sourceColumn, Index targetColumnIndex)
{
  if constexpr (std::is_integral_v<Entry_range_or_column_index>) {
    _get_column(targetColumnIndex) += _get_column(sourceColumn);
  } else {
    _get_column(targetColumnIndex) += sourceColumn;
  }
}

template <class Master_matrix>
template <class Entry_range_or_column_index>
inline void Base_matrix<Master_matrix>::multiply_target_and_add_to(const Entry_range_or_column_index& sourceColumn,
                                                                   const Field_element& coefficient,
                                                                   Index targetColumnIndex)
{
  if constexpr (std::is_integral_v<Entry_range_or_column_index>) {
    _get_column(targetColumnIndex).multiply_target_and_add(coefficient, _get_column(sourceColumn));
  } else {
    _get_column(targetColumnIndex).multiply_target_and_add(coefficient, sourceColumn);
  }
}

template <class Master_matrix>
template <class Entry_range_or_column_index>
inline void Base_matrix<Master_matrix>::multiply_source_and_add_to(const Field_element& coefficient,
                                                                   const Entry_range_or_column_index& sourceColumn,
                                                                   Index targetColumnIndex)
{
  if constexpr (std::is_integral_v<Entry_range_or_column_index>) {
    _get_column(targetColumnIndex).multiply_source_and_add(_get_column(sourceColumn), coefficient);
  } else {
    _get_column(targetColumnIndex).multiply_source_and_add(sourceColumn, coefficient);
  }
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_entry(Index columnIndex, Index rowIndex)
{
  _get_column(columnIndex).clear(_get_real_row_index(rowIndex));
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_column(Index columnIndex)
{
  _get_column(columnIndex).clear();
}

template <class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_entry(Index columnIndex, Index rowIndex) const
{
  return !(_get_column(columnIndex).is_non_zero(_get_real_row_index(rowIndex)));
}

template <class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_column(Index columnIndex)
{
  return _get_column(columnIndex).is_empty();
}

template <class Master_matrix>
inline Base_matrix<Master_matrix>& Base_matrix<Master_matrix>::operator=(const Base_matrix& other)
{
  if (this == &other) return *this;

  Swap_opt::operator=(other);
  RA_opt::operator=(other);
  matrix_.clear();
  nextInsertIndex_ = other.nextInsertIndex_;
  colSettings_ = other.colSettings_;

  matrix_.reserve(other.matrix_.size());
  for (const auto& cont : other.matrix_) {
    _container_insert(cont);
  }

  return *this;
}

template <class Master_matrix>
inline Base_matrix<Master_matrix>& Base_matrix<Master_matrix>::operator=(Base_matrix&& other) noexcept
{
  if (this == &other) return *this;

  Swap_opt::operator=(std::move(other));
  RA_opt::operator=(std::move(other));

  matrix_ = std::move(other.matrix_);
  nextInsertIndex_ = std::exchange(other.nextInsertIndex_, 0);
  colSettings_ = std::exchange(other.colSettings_, nullptr);

  return *this;
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::print()
{
  _orderRowsIfNecessary();
  std::cout << "Base_matrix:\n";
  for (Index i = 0; i < nextInsertIndex_; ++i) {
    const Column& col = matrix_[i];
    for (const auto& e : col.get_content(nextInsertIndex_)) {
      if (e == 0U)
        std::cout << "- ";
      else
        std::cout << e << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  if constexpr (Master_matrix::Option_list::has_row_access) {
    std::cout << "Row Matrix:\n";
    for (Index i = 0; i < nextInsertIndex_; ++i) {
      const auto& row = RA_opt::get_row(i);
      for (const auto& entry : row) {
        std::cout << entry.get_column_index() << " ";
      }
      std::cout << "(" << i << ")\n";
    }
    std::cout << "\n";
  }
}

template <class Master_matrix>
template <class Container, class>
inline void Base_matrix<Master_matrix>::_insert(const Container& column, Index columnIndex, Dimension dim)
{
  _orderRowsIfNecessary();

  // resize of containers when necessary:
  Index pivot = 0;
  if (column.begin() != column.end()) {
    // first, compute pivot of `column`
    pivot = Master_matrix::get_row_index(*std::prev(column.end()));
    // row container
    if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_removable_rows)
      RA_opt::_resize(pivot);
  }

  // row swap map containers
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
      for (const auto& id : column) {
        Swap_opt::_initialize_row_index(Master_matrix::get_row_index(id));
      }
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
      Swap_opt::_initialize_row_index(pivot);
    }
    // column container
    if constexpr (!Master_matrix::Option_list::has_row_access) {
      if (matrix_.size() <= columnIndex) {
        matrix_.resize(columnIndex + 1);
      }
    }
  }

  _container_insert(column, columnIndex, dim);
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::_insert(Index idx, Field_element e, Index columnIndex, Dimension dim)
{
  _orderRowsIfNecessary();

  // resize of containers when necessary:
  if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_removable_rows)
    RA_opt::_resize(idx);

  // row swap map containers
  if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
    Swap_opt::_initialize_row_index(idx);
  }

  // column container
  if constexpr (!Master_matrix::Option_list::has_map_column_container && !Master_matrix::Option_list::has_row_access) {
    if (matrix_.size() <= columnIndex) {
      matrix_.resize(columnIndex + 1);
    }
  }

  _container_insert(idx, e, columnIndex, dim);
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::_orderRowsIfNecessary()
{
  if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
    if (Swap_opt::_row_were_swapped()) Swap_opt::_orderRows();
  }
}

template <class Master_matrix>
inline const typename Base_matrix<Master_matrix>::Column& Base_matrix<Master_matrix>::_get_column(
    Index columnIndex) const
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.at(columnIndex);
  } else {
    return matrix_[columnIndex];
  }
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::Column& Base_matrix<Master_matrix>::_get_column(Index columnIndex)
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.at(columnIndex);
  } else {
    return matrix_[columnIndex];
  }
}

template <class Master_matrix>
inline typename Base_matrix<Master_matrix>::Index Base_matrix<Master_matrix>::_get_real_row_index(Index rowIndex) const
{
  if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
    return Swap_opt::_get_row_index(rowIndex);
  } else {
    return rowIndex;
  }
}

template <class Master_matrix>
template <class Container>
inline void Base_matrix<Master_matrix>::_container_insert(const Container& column, Index pos, Dimension dim)
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.try_emplace(pos, Column(pos, column, dim, RA_opt::_get_rows_ptr(), colSettings_));
    } else {
      matrix_.try_emplace(pos, Column(column, dim, colSettings_));
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.emplace_back(pos, column, dim, RA_opt::_get_rows_ptr(), colSettings_);
    } else {
      matrix_[pos] = Column(column, dim, colSettings_);
    }
  }
}

template <class Master_matrix>
inline void Base_matrix<Master_matrix>::_container_insert(Index idx, [[maybe_unused]] Field_element e, Index pos, Dimension dim){
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      if constexpr (Master_matrix::Option_list::is_z2){
        matrix_.try_emplace(pos, Column(pos, idx, dim, RA_opt::_get_rows_ptr(), colSettings_));
      } else {
        matrix_.try_emplace(pos, Column(pos, idx, e, dim, RA_opt::_get_rows_ptr(), colSettings_));
      }
    } else {
      if constexpr (Master_matrix::Option_list::is_z2){
        matrix_.try_emplace(pos, Column(idx, dim, colSettings_));
      } else {
        matrix_.try_emplace(pos, Column(idx, e, dim, colSettings_));
      }
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      if constexpr (Master_matrix::Option_list::is_z2){
        matrix_.emplace_back(pos, idx, dim, RA_opt::_get_rows_ptr(), colSettings_);
      } else {
        matrix_.emplace_back(pos, idx, e, dim, RA_opt::_get_rows_ptr(), colSettings_);
      }
    } else {
      if constexpr (Master_matrix::Option_list::is_z2){
        matrix_[pos] = Column(idx, dim, colSettings_);
      } else {
        matrix_[pos] = Column(idx, e, dim, colSettings_);
      }
    }
  }
}

template <class Master_matrix>
template <class ColumnIterator> // Pair (pos,Column) if has_map_column_container, Column otherwise
inline void Base_matrix<Master_matrix>::_container_insert(const ColumnIterator& rep)
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    const auto& col = rep.second;
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.try_emplace(rep.first, Column(col, col.get_column_index(), RA_opt::_get_rows_ptr(), colSettings_));
    } else {
      matrix_.try_emplace(rep.first, Column(col, colSettings_));
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.emplace_back(rep, rep.get_column_index(), RA_opt::_get_rows_ptr(), colSettings_);
    } else {
      matrix_.emplace_back(rep, colSettings_);
    }
  }
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_BASE_MATRIX_H
