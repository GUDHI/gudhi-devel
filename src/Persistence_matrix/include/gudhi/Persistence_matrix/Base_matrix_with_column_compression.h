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
 * @file Base_matrix_with_column_compression.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Base_matrix_with_column_compression class.
 */

#ifndef PM_BASE_MATRIX_COMPRESSION_H
#define PM_BASE_MATRIX_COMPRESSION_H

#include <iostream>   //print() only
#include <vector>
#include <utility>    //std::swap, std::move & std::exchange

#include <boost/intrusive/set.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include <gudhi/Simple_object_pool.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Base_matrix_with_column_compression Base_matrix_with_column_compression.h \
 * gudhi/Persistence_matrix/Base_matrix_with_column_compression.h
 * @ingroup persistence_matrix
 *
 * @brief A @ref basematrix "base matrix" (also see @ref Base_matrix), but with column compression. That is, all
 * identical columns in the matrix are compressed together as the same column. For matrices with a lot of redundant
 * columns, this will save a lot of space. Also, any addition made onto a column will be performed at the same time
 * on all other identical columns, which is an advantage for the cohomology algorithm for example.
 * 
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Base_matrix_with_column_compression : protected Master_matrix::Matrix_row_access_option 
{
 public:
  using Index = typename Master_matrix::Index;                          /**< Container index type. */
  using Dimension = typename Master_matrix::Dimension;                  /**< Dimension value type. */
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element = typename Master_matrix::Element;                /**< Field element value type. */
  using Row = typename Master_matrix::Row;                              /**< Row type,
                                                                             only necessary with row access option. */
  using Entry_constructor = typename Master_matrix::Entry_constructor;  /**< Factory of @ref Entry classes. */
  using Column_settings = typename Master_matrix::Column_settings;      /**< Structure giving access to the columns to
                                                                             necessary external classes. */

  /**
   * @brief Type for columns. Only one for each "column class" is explicitly constructed.
   */
  class Column
      : public Master_matrix::Column,
        public boost::intrusive::set_base_hook<boost::intrusive::link_mode<boost::intrusive::normal_link> > 
  {
   public:
    using Base = typename Master_matrix::Column;

    Column(Column_settings* colSettings = nullptr)
        : Base(colSettings) {}
    template <class Container>
    Column(const Container& nonZeroRowIndices, Column_settings* colSettings)
        : Base(nonZeroRowIndices, colSettings) {}
    template <class Container, class Row_container>
    Column(Index columnIndex, const Container& nonZeroRowIndices, Row_container& rowContainer,
                Column_settings* colSettings)
        : Base(columnIndex, nonZeroRowIndices, rowContainer, colSettings) {}
    template <class Container>
    Column(const Container& nonZeroRowIndices, Dimension dimension, Column_settings* colSettings)
        : Base(nonZeroRowIndices, dimension, colSettings) {}
    template <class Container, class Row_container>
    Column(Index columnIndex, const Container& nonZeroRowIndices, Dimension dimension,
                Row_container& rowContainer, Column_settings* colSettings)
        : Base(columnIndex, nonZeroRowIndices, dimension, rowContainer, colSettings) {}
    Column(const Column& column, Column_settings* colSettings = nullptr)
        : Base(static_cast<const Base&>(column), colSettings) {}
    template <class Row_container>
    Column(const Column& column, Index columnIndex, Row_container& rowContainer,
                Column_settings* colSettings = nullptr)
        : Base(static_cast<const Base&>(column), columnIndex, rowContainer, colSettings) {}
    Column(Column&& column) noexcept : Base(std::move(static_cast<Base&>(column))) {}

    //TODO: is it possible to make this work?
    // template <class... U>
    // Column(U&&... u) : Base(std::forward<U>(u)...) {}

    Index get_rep() const { return rep_; }
    void set_rep(const Index& rep) { rep_ = rep; }

    struct Hasher {
      size_t operator()(const Column& c) const { return std::hash<Base>()(static_cast<Base>(c)); }
    };

   private:
    Index rep_; /**< Index in the union-find of the root of the set representing this column class. */
  };

  /**
   * @brief Constructs an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Base_matrix_with_column_compression(Column_settings* colSettings);
  /**
   * @brief Constructs a matrix from the given ordered columns. The columns are inserted in the given order.
   * If no identical column already existed, a copy of the column is stored. If an identical one existed, no new
   * column is constructed and the relationship between the two is registered in an union-find structure.
   * 
   * @tparam Container Range type for @ref Matrix::Entry_representative ranges.
   * Assumed to have a begin(), end() and size() method.
   * @param columns A vector of @ref Matrix::Entry_representative ranges to construct the columns from.
   * The content of the ranges are assumed to be sorted by increasing ID value.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  template <class Container>
  Base_matrix_with_column_compression(const std::vector<Container>& columns, 
                                      Column_settings* colSettings);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Base_matrix_with_column_compression(unsigned int numberOfColumns, 
                                      Column_settings* colSettings);
  /**
   * @brief Copy constructor. If @p colSettings is not a null pointer, its value is kept
   * instead of the one in the copied matrix.
   * 
   * @param matrixToCopy Matrix to copy.
   * @param colSettings Either a pointer to an existing setting structure for the columns or a null pointer.
   * The structure should contain all the necessary external classes specifically necessary for the choosen column type,
   * such as custom allocators. If null pointer, the pointer stored in @p matrixToCopy is used instead.
   */
  Base_matrix_with_column_compression(const Base_matrix_with_column_compression& matrixToCopy,
                                      Column_settings* colSettings = nullptr);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Base_matrix_with_column_compression(Base_matrix_with_column_compression&& other) noexcept;
  /**
   * @brief Destructor.
   */
  ~Base_matrix_with_column_compression();

  /**
   * @brief Inserts a new ordered column at the end of the matrix by copying the given range of
   * @ref Matrix::Entry_representative. The content of the range is assumed to be sorted by increasing ID value. 
   * 
   * @tparam Container Range of @ref Matrix::Entry_representative. Assumed to have a begin(), end() and size() method.
   * @param column Range of @ref Matrix::Entry_representative from which the column has to be constructed. Assumed to be
   * ordered by increasing ID value. 
   */
  template <class Container>
  void insert_column(const Container& column);
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
  void insert_boundary(const Boundary_range& boundary, Dimension dim = -1);
  /**
   * @brief Returns the column at the given @ref MatIdx index.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   *
   * Remark: the method it-self is not const, because of the path compression optimization of the union-find structure,
   * when a column is looked up. 
   * 
   * @param columnIndex @ref MatIdx index of the column to return.
   * @return Const reference to the column.
   */
  const Column& get_column(Index columnIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access is true.
   * Returns the row at the given @ref rowindex "row index" of the compressed matrix.
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * Note that the row will be from the compressed matrix, that is, the one with only unique columns.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return.
   * @return Const reference to the row.
   */
  const Row& get_row(Index rowIndex) const;
  /**
   * @brief If @ref PersistenceMatrixOptions::has_row_access and @ref PersistenceMatrixOptions::has_removable_rows
   * are true: assumes that the row is empty and removes it. Otherwise, does nothing.
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
   * @brief Returns the current number of columns in the matrix, counting also the redundant columns.
   * 
   * @return The number of columns.
   */
  Index get_number_of_columns() const;

  /**
   * @brief Adds column represented by @p sourceColumn onto the column at @p targetColumnIndex in the matrix.
   *
   * The representatives of redundant columns are summed together, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
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
   * The representatives of redundant columns are summed together, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
   * 
   * @tparam Entry_range_or_column_index Either a range of @ref Entry with a begin() and end() method,
   * or any integer type.
   * @param sourceColumn Either a @ref Entry range or the @ref MatIdx index of the column to add.
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
   * The representatives of redundant columns are summed together, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
   * 
   * @tparam Entry_range_or_column_index Either a range of @ref Entry with a begin() and end() method,
   * or any integer type.
   * @param coefficient Value to multiply.
   * @param sourceColumn Either a @ref Entry range or the @ref MatIdx index of the column to add.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Entry_range_or_column_index>
  void multiply_source_and_add_to(const Field_element& coefficient, 
                                  const Entry_range_or_column_index& sourceColumn,
                                  Index targetColumnIndex);

  /**
   * @brief Indicates if the entry at given coordinates has value zero.
   * 
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @return true If the entry has value zero.
   * @return false Otherwise.
   */
  bool is_zero_entry(Index columnIndex, Index rowIndex);
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
   * @param colSettings Pointer to the entry factory.
   */
  void reset(Column_settings* colSettings) {
    columnToRep_.clear_and_dispose(Delete_disposer(this));
    columnClasses_ = boost::disjoint_sets_with_storage<>();
    repToColumn_.clear();
    nextColumnIndex_ = 0;
    colSettings_ = colSettings;
  }

  /**
   * @brief Assign operator.
   */
  Base_matrix_with_column_compression& operator=(const Base_matrix_with_column_compression& other);
  /**
   * @brief Swap operator.
   */
  friend void swap(Base_matrix_with_column_compression& matrix1, Base_matrix_with_column_compression& matrix2) {
    matrix1.columnToRep_.swap(matrix2.columnToRep_);
    std::swap(matrix1.columnClasses_, matrix2.columnClasses_);
    matrix1.repToColumn_.swap(matrix2.repToColumn_);
    std::swap(matrix1.nextColumnIndex_, matrix2.nextColumnIndex_);
    std::swap(matrix1.colSettings_, matrix2.colSettings_);
    std::swap(matrix1.columnPool_, matrix2.columnPool_);

    if constexpr (Master_matrix::Option_list::has_row_access) {
      swap(static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix1),
           static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix2));
    }
  }

  void print();  // for debug

 private:
  /**
   * @brief The disposer object function for boost intrusive container
   */
  struct Delete_disposer {
    Delete_disposer(Base_matrix_with_column_compression* matrix) : matrix_(matrix) {}

    void operator()(Column* delete_this) { matrix_->columnPool_->destroy(delete_this); }

    Base_matrix_with_column_compression* matrix_;
  };

  using RA_opt = typename Master_matrix::Matrix_row_access_option;
  using Col_dict = boost::intrusive::set<Column, boost::intrusive::constant_time_size<false> >;

  Col_dict columnToRep_;                              /**< Map from a column to the index of its representative. */
  boost::disjoint_sets_with_storage<> columnClasses_; /**< Union-find structure,
                                                           where two columns in the same set are identical. */
  std::vector<Column*> repToColumn_;                  /**< Map from the representative index to
                                                           the representative Column. */
  Index nextColumnIndex_;                             /**< Next unused column index. */
  Column_settings* colSettings_;                      /**< Entry factory. */
  /**
   * @brief Column factory. Has to be a pointer as Simple_object_pool is not swappable, so their addresses have to be
   * exchanged instead.
   */
  std::unique_ptr<Simple_object_pool<Column> > columnPool_;
  inline static const Column empty_column_;           /**< Representative for empty columns. */

  void _insert_column(Index columnIndex);
  void _insert_double_column(Index columnIndex, typename Col_dict::iterator& doubleIt);
};

template <class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(
    Column_settings* colSettings)
    : RA_opt(), nextColumnIndex_(0), colSettings_(colSettings), columnPool_(new Simple_object_pool<Column>)
{}

template <class Master_matrix>
template <class Container>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(
    const std::vector<Container>& columns, Column_settings* colSettings)
    : RA_opt(columns.size()),
      columnClasses_(columns.size()),
      repToColumn_(columns.size(), nullptr),
      nextColumnIndex_(0),
      colSettings_(colSettings),
      columnPool_(new Simple_object_pool<Column>)
{
  for (const Container& c : columns) {
    insert_column(c);
  }
}

template <class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(
    unsigned int numberOfColumns, Column_settings* colSettings)
    : RA_opt(numberOfColumns),
      columnClasses_(numberOfColumns),
      repToColumn_(numberOfColumns, nullptr),
      nextColumnIndex_(0),
      colSettings_(colSettings),
      columnPool_(new Simple_object_pool<Column>)
{}

template <class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(
    const Base_matrix_with_column_compression& matrixToCopy, Column_settings* colSettings)
    : RA_opt(static_cast<const RA_opt&>(matrixToCopy)),
      columnClasses_(matrixToCopy.columnClasses_),
      repToColumn_(matrixToCopy.repToColumn_.size(), nullptr),
      nextColumnIndex_(0),
      colSettings_(colSettings == nullptr ? matrixToCopy.colSettings_ : colSettings),
      columnPool_(new Simple_object_pool<Column>)
{
  for (const Column* col : matrixToCopy.repToColumn_) {
    if (col != nullptr) {
      if constexpr (Master_matrix::Option_list::has_row_access) {
        repToColumn_[nextColumnIndex_] =
            columnPool_->construct(*col, col->get_column_index(), RA_opt::rows_, colSettings_);
      } else {
        repToColumn_[nextColumnIndex_] = columnPool_->construct(*col, colSettings_);
      }
      columnToRep_.insert(columnToRep_.end(), *repToColumn_[nextColumnIndex_]);
      repToColumn_[nextColumnIndex_]->set_rep(nextColumnIndex_);
    }
    nextColumnIndex_++;
  }
}

template <class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(
    Base_matrix_with_column_compression&& other) noexcept
    : RA_opt(std::move(static_cast<RA_opt&>(other))),
      columnToRep_(std::move(other.columnToRep_)),
      columnClasses_(std::move(other.columnClasses_)),
      repToColumn_(std::move(other.repToColumn_)),
      nextColumnIndex_(std::exchange(other.nextColumnIndex_, 0)),
      colSettings_(std::exchange(other.colSettings_, nullptr)),
      columnPool_(std::exchange(other.columnPool_, nullptr))
{}

template <class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::~Base_matrix_with_column_compression() 
{
  columnToRep_.clear_and_dispose(Delete_disposer(this));
}

template <class Master_matrix>
template <class Container>
inline void Base_matrix_with_column_compression<Master_matrix>::insert_column(const Container& column) 
{
  insert_boundary(column);
}

template <class Master_matrix>
template <class Boundary_range>
inline void Base_matrix_with_column_compression<Master_matrix>::insert_boundary(const Boundary_range& boundary,
                                                                                Dimension dim) 
{
  // handles a dimension which is not actually stored.
  // TODO: verify if this is not a problem for the cohomology algorithm and if yes,
  // change how Column_dimension_option is choosen.
  if (dim == -1) dim = boundary.size() == 0 ? 0 : boundary.size() - 1;

  if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_removable_rows) {
    if (boundary.begin() != boundary.end()) {
      Index pivot;
      if constexpr (Master_matrix::Option_list::is_z2) {
        pivot = *std::prev(boundary.end());
      } else {
        pivot = std::prev(boundary.end())->first;
      }
      if (RA_opt::rows_->size() <= pivot) RA_opt::rows_->resize(pivot + 1);
    }
  }

  if (repToColumn_.size() == nextColumnIndex_) {
    // could perhaps be avoided, if find_set returns something special when it does not find
    columnClasses_.link(nextColumnIndex_, nextColumnIndex_);
    if constexpr (Master_matrix::Option_list::has_row_access) {
      repToColumn_.push_back(
          columnPool_->construct(nextColumnIndex_, boundary, dim, RA_opt::rows_, colSettings_));
    } else {
      repToColumn_.push_back(columnPool_->construct(boundary, dim, colSettings_));
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      repToColumn_[nextColumnIndex_] =
          columnPool_->construct(nextColumnIndex_, boundary, dim, RA_opt::rows_, colSettings_);
    } else {
      repToColumn_[nextColumnIndex_] = columnPool_->construct(boundary, dim, colSettings_);
    }
  }
  _insert_column(nextColumnIndex_);

  nextColumnIndex_++;
}

template <class Master_matrix>
inline const typename Base_matrix_with_column_compression<Master_matrix>::Column&
Base_matrix_with_column_compression<Master_matrix>::get_column(Index columnIndex) 
{
  auto col = repToColumn_[columnClasses_.find_set(columnIndex)];
  if (col == nullptr) return empty_column_;
  return *col;
}

template <class Master_matrix>
inline const typename Base_matrix_with_column_compression<Master_matrix>::Row&
Base_matrix_with_column_compression<Master_matrix>::get_row(Index rowIndex) const 
{
  static_assert(Master_matrix::Option_list::has_row_access, "Row access has to be enabled for this method.");

  return RA_opt::get_row(rowIndex);
}

template <class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::erase_empty_row(Index rowIndex) 
{
  if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows) {
    RA_opt::erase_empty_row(rowIndex);
  }
}

template <class Master_matrix>
inline typename Base_matrix_with_column_compression<Master_matrix>::Index
Base_matrix_with_column_compression<Master_matrix>::get_number_of_columns() const 
{
  return nextColumnIndex_;
}

template <class Master_matrix>
template <class Entry_range_or_column_index>
inline void Base_matrix_with_column_compression<Master_matrix>::add_to(const Entry_range_or_column_index& sourceColumn,
                                                                       Index targetColumnIndex) 
{
  // handle case where targetRep == sourceRep?
  Index targetRep = columnClasses_.find_set(targetColumnIndex);
  Column& target = *repToColumn_[targetRep];
  columnToRep_.erase(target);
  if constexpr (std::is_integral_v<Entry_range_or_column_index>) {
    target += get_column(sourceColumn);
  } else {
    target += sourceColumn;
  }
  _insert_column(targetRep);
}

template <class Master_matrix>
template <class Entry_range_or_column_index>
inline void Base_matrix_with_column_compression<Master_matrix>::multiply_target_and_add_to(
    const Entry_range_or_column_index& sourceColumn, const Field_element& coefficient, Index targetColumnIndex) 
{
  // handle case where targetRep == sourceRep?
  Index targetRep = columnClasses_.find_set(targetColumnIndex);
  Column& target = *repToColumn_[targetRep];
  columnToRep_.erase(target);
  if constexpr (std::is_integral_v<Entry_range_or_column_index>) {
    target.multiply_target_and_add(coefficient, get_column(sourceColumn));
  } else {
    target.multiply_target_and_add(coefficient, sourceColumn);
  }
  _insert_column(targetRep);
}

template <class Master_matrix>
template <class Entry_range_or_column_index>
inline void Base_matrix_with_column_compression<Master_matrix>::multiply_source_and_add_to(
    const Field_element& coefficient, const Entry_range_or_column_index& sourceColumn, Index targetColumnIndex) 
{
  // handle case where targetRep == sourceRep?
  Index targetRep = columnClasses_.find_set(targetColumnIndex);
  Column& target = *repToColumn_[targetRep];
  columnToRep_.erase(target);
  if constexpr (std::is_integral_v<Entry_range_or_column_index>) {
    target.multiply_source_and_add(get_column(sourceColumn), coefficient);
  } else {
    target.multiply_source_and_add(sourceColumn, coefficient);
  }
  _insert_column(targetRep);
}

template <class Master_matrix>
inline bool Base_matrix_with_column_compression<Master_matrix>::is_zero_entry(Index columnIndex, Index rowIndex) 
{
  auto col = repToColumn_[columnClasses_.find_set(columnIndex)];
  if (col == nullptr) return true;
  return !col->is_non_zero(rowIndex);
}

template <class Master_matrix>
inline bool Base_matrix_with_column_compression<Master_matrix>::is_zero_column(Index columnIndex) 
{
  auto col = repToColumn_[columnClasses_.find_set(columnIndex)];
  if (col == nullptr) return true;
  return col->is_empty();
}

template <class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>&
Base_matrix_with_column_compression<Master_matrix>::operator=(const Base_matrix_with_column_compression& other) 
{
  for (auto col : repToColumn_) {
    if (col != nullptr) {
      columnPool_->destroy(col);
      col = nullptr;
    }
  }
  RA_opt::operator=(other);
  columnClasses_ = other.columnClasses_;
  columnToRep_.reserve(other.columnToRep_.size());
  repToColumn_.resize(other.repToColumn_.size(), nullptr);
  nextColumnIndex_ = 0;
  colSettings_ = other.colSettings_;
  for (const Column* col : other.repToColumn_) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      repToColumn_[nextColumnIndex_] =
          columnPool_->construct(*col, col->get_column_index(), RA_opt::rows_, colSettings_);
    } else {
      repToColumn_[nextColumnIndex_] = columnPool_->construct(*col, colSettings_);
    }
    columnToRep_.insert(columnToRep_.end(), *repToColumn_[nextColumnIndex_]);
    repToColumn_[nextColumnIndex_]->set_rep(nextColumnIndex_);
    nextColumnIndex_++;
  }
  return *this;
}

template <class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::print() 
{
  std::cout << "Compressed_matrix:\n";
  for (Column& col : columnToRep_) {
    for (auto e : col->get_content(nextColumnIndex_)) {
      if (e == 0u)
        std::cout << "- ";
      else
        std::cout << e << " ";
    }
    std::cout << "(";
    for (Index i = 0; i < nextColumnIndex_; ++i) {
      if (columnClasses_.find_set(i) == col.get_rep()) std::cout << i << " ";
    }
    std::cout << ")\n";
  }
  std::cout << "\n";
  std::cout << "Row Matrix:\n";
  for (Index i = 0; i < RA_opt::rows_->size(); ++i) {
    const Row& row = RA_opt::rows_[i];
    for (const auto& entry : row) {
      std::cout << entry.get_column_index() << " ";
    }
    std::cout << "(" << i << ")\n";
  }
  std::cout << "\n";
}

template <class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::_insert_column(Index columnIndex) 
{
  Column& col = *repToColumn_[columnIndex];

  if (col.is_empty()) {
    columnPool_->destroy(&col);  // delete curr_col;
    repToColumn_[columnIndex] = nullptr;
    return;
  }

  col.set_rep(columnIndex);
  auto res = columnToRep_.insert(col);
  if (res.first->get_rep() != columnIndex) {  //if true, then redundant column
    _insert_double_column(columnIndex, res.first);
  }
}

template <class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::_insert_double_column(
    Index columnIndex, typename Col_dict::iterator& doubleIt) 
{
  Index doubleRep = doubleIt->get_rep();
  columnClasses_.link(columnIndex, doubleRep);  // both should be representatives
  Index newRep = columnClasses_.find_set(columnIndex);

  columnPool_->destroy(repToColumn_[columnIndex]);
  repToColumn_[columnIndex] = nullptr;

  if (newRep == columnIndex) {
    std::swap(repToColumn_[doubleRep], repToColumn_[columnIndex]);
    doubleIt->set_rep(columnIndex);
  }
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_BASE_MATRIX_COMPRESSION_H
