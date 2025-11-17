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
 * @file Boundary_matrix.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Boundary_matrix class.
 */

#ifndef PM_BOUNDARY_MATRIX_H
#define PM_BOUNDARY_MATRIX_H

#include <cassert>
#include <iostream>  //print() only
#include <vector>
#include <utility>  //std::swap, std::move & std::exchange

namespace Gudhi {
namespace persistence_matrix {

// TODO: factorize/inherit/compose with Base_matrix?
/**
 * @class Boundary_matrix Boundary_matrix.h gudhi/Persistence_matrix/Boundary_matrix.h
 * @ingroup persistence_matrix
 *
 * @brief %Matrix structure to store the ordered @ref boundarymatrix "boundary matrix" \f$ R \f$ of a filtered complex
 * in order to compute its persistent homology. Provides an access to its columns and rows as well as the possibility
 * to remove the last cells of the filtration while maintaining a valid barcode.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Boundary_matrix : public Master_matrix::Matrix_dimension_option,
                        public Master_matrix::template Base_swap_option<Boundary_matrix<Master_matrix> >,
                        public Master_matrix::Base_pairing_option,
                        protected Master_matrix::Matrix_row_access_option
{
 private:
  using Dim_opt = typename Master_matrix::Matrix_dimension_option;
  using Swap_opt = typename Master_matrix::template Base_swap_option<Boundary_matrix<Master_matrix> >;
  using Pair_opt = typename Master_matrix::Base_pairing_option;
  using RA_opt = typename Master_matrix::Matrix_row_access_option;

  static constexpr bool activeDimOption_ =
      Master_matrix::Option_list::has_matrix_maximal_dimension_access || Master_matrix::maxDimensionIsNeeded;
  static constexpr bool activeSwapOption_ =
      Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update;
  static constexpr bool activePairingOption_ = Master_matrix::Option_list::has_column_pairings &&
                                               !Master_matrix::Option_list::has_vine_update &&
                                               !Master_matrix::Option_list::can_retrieve_representative_cycles;

 public:
  using Index = typename Master_matrix::Index;         /**< Container index type. */
  using ID_index = typename Master_matrix::ID_index;   /**< @ref IDIdx index type. */
  using Dimension = typename Master_matrix::Dimension; /**< Dimension value type. */
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element = typename Master_matrix::Element;               /**< Type of an field element. */
  using Column = typename Master_matrix::Column;                       /**< Column type. */
  using Boundary = typename Master_matrix::Boundary;                   /**< Type of an input column. */
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
  Boundary_matrix(Column_settings* colSettings);
  /**
   * @brief Constructs a new matrix from the given ranges of @ref Matrix::Entry_representative. Each range corresponds
   * to a column  (the order of the ranges are preserved). The content of the ranges is assumed to be sorted by
   * increasing IDs. The IDs of the simplices are also assumed to be consecutive, ordered by filtration value, starting
   * with 0.
   *
   * @tparam Boundary_range Range type for @ref Matrix::Entry_representative ranges.
   * Assumed to have a begin(), end() and size() method.
   * @param orderedBoundaries Range of boundaries: @p orderedBoundaries is interpreted as a boundary matrix of a
   * filtered **simplicial** complex, whose boundaries are ordered by filtration order.
   * Therefore, `orderedBoundaries[i]` should store the boundary of the \f$ i^{th} \f$ simplex in the filtration,
   * as an ordered list of indices of its facets (again those indices correspond to their respective position
   * in the matrix). That is why the indices of the simplices are assumed to be consecutive and starting with 0
   * (an empty boundary is interpreted as a vertex boundary and not as a non existing simplex).
   * All dimensions up to the maximal dimension of interest have to be present. If only a higher dimension is of
   * interest and not everything should be stored, then use the @ref insert_boundary method instead
   * (after creating the matrix with the
   * @ref Boundary_matrix(unsigned int numberOfColumns, Column_settings* colSettings)
   * constructor preferably).
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  template <class Boundary_range = Boundary>
  Boundary_matrix(const std::vector<Boundary_range>& orderedBoundaries, Column_settings* colSettings);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   *
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  Boundary_matrix(unsigned int numberOfColumns, Column_settings* colSettings);
  /**
   * @brief Copy constructor. If @p colSettings is not a null pointer, its value is kept
   * instead of the one in the copied matrix.
   *
   * @param matrixToCopy Matrix to copy.
   * @param colSettings Either a pointer to an existing setting structure for the columns or a null pointer.
   * The structure should contain all the necessary external classes specifically necessary for the chosen column type,
   * such as custom allocators. If null pointer, the pointer stored in @p matrixToCopy is used instead.
   */
  Boundary_matrix(const Boundary_matrix& matrixToCopy, Column_settings* colSettings = nullptr);
  /**
   * @brief Move constructor.
   *
   * @param other Matrix to move.
   */
  Boundary_matrix(Boundary_matrix&& other) noexcept;

  ~Boundary_matrix() = default;

  /**
   * @brief Inserts at the end of the matrix a new ordered column corresponding to the given boundary.
   * This means that it is assumed that this method is called on boundaries in the order of the filtration.
   * It also assumes that the cells in the given boundary are identified by their relative position in the filtration,
   * starting at 0. If it is not the case, use the other
   * @ref insert_boundary(ID_index cellIndex, const Boundary_range& boundary, Dimension dim) "insert_boundary"
   * instead by indicating the cell ID used in the boundaries when the cell is inserted.
   *
   * Different to the constructor, the boundaries do not have to come from a simplicial complex, but also from
   * a more general entry complex. This includes cubical complexes or Morse complexes for example.
   *
   * At the insertion, the boundary will be copied as is. The column will only be reduced later when the barcode
   * is requested in order to apply some optimizations with the additional knowledge. Hence, the barcode will also
   * not be updated, so call @ref Base_pairing::get_current_barcode "get_current_barcode" only when the matrix is
   * complete.
   *
   * @tparam Boundary_range Range of @ref Matrix::Entry_representative. Assumed to have a begin(), end() and size()
   * method.
   * @param boundary Boundary generating the new column. The content should be ordered by ID.
   * @param dim Dimension of the cell whose boundary is given. If the complex is simplicial,
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return The @ref MatIdx index of the inserted boundary.
   */
  template <class Boundary_range = Boundary>
  Index insert_boundary(const Boundary_range& boundary,
                        Dimension dim = Master_matrix::template get_null_value<Dimension>());
  /**
   * @brief It does the same as the other version, but allows the boundary cells to be identified without restrictions
   * except that all IDs have to be strictly increasing in the order of filtration. Note that you should avoid then
   * to use the other insertion method to avoid overwriting IDs.
   *
   * As a cell has to be inserted before one of its cofaces in a valid filtration (recall that it is assumed that
   * the cells are inserted by order of filtration), it is sufficient to indicate the ID of the cell being inserted.
   *
   * @tparam Boundary_range Range of @ref Matrix::Entry_representative. Assumed to have a begin(), end() and size()
   * method.
   * @param cellIndex @ref IDIdx index to use to identify the new cell.
   * @param boundary Boundary generating the new column. The indices of the boundary have to correspond to the
   * @p cellIndex values of precedent calls of the method for the corresponding cells and should be ordered in
   * increasing order.
   * @param dim Dimension of the cell whose boundary is given. If the complex is simplicial,
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return The @ref MatIdx index of the inserted boundary.
   */
  template <class Boundary_range = Boundary>
  Index insert_boundary(ID_index cellIndex,
                        const Boundary_range& boundary,
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
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns is true.
   * Removes the last cell in the filtration from the matrix and updates the barcode if this one was already computed.
   *
   * @return The pivot of the removed cell.
   */
  Index remove_last();
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
   * @brief Returns the dimension of the given column.
   *
   * @param columnIndex @ref MatIdx index of the column representing the cell.
   * @return Dimension of the cell.
   */
  Dimension get_column_dimension(Index columnIndex) const;

  /**
   * @brief Adds column at @p sourceColumnIndex onto the column at @p targetColumnIndex in the matrix.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of a
   * boundary matrix of a filtered complex. For example, a right-to-left addition could corrupt the computation
   * of the barcode if done blindly. So should be used with care.
   *
   * @param sourceColumnIndex @ref MatIdx index of the source column.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  void add_to(Index sourceColumnIndex, Index targetColumnIndex);
  /**
   * @brief Multiplies the target column with the coefficient and then adds the source column to it.
   * That is: `targetColumn = (targetColumn * coefficient) + sourceColumn`.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of a
   * boundary matrix of a filtered complex. For example, a right-to-left addition could corrupt the computation
   * of the barcode if done blindly. So should be used with care.
   *
   * @param sourceColumnIndex @ref MatIdx index of the source column.
   * @param coefficient Value to multiply.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  void multiply_target_and_add_to(Index sourceColumnIndex, const Field_element& coefficient, Index targetColumnIndex);
  /**
   * @brief Multiplies the source column with the coefficient before adding it to the target column.
   * That is: `targetColumn += (coefficient * sourceColumn)`. The source column will **not** be modified.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of a
   * boundary matrix of a filtered complex. For example, a right-to-left addition could corrupt the computation
   * of the barcode if done blindly. So should be used with care.
   *
   * @param coefficient Value to multiply.
   * @param sourceColumnIndex @ref MatIdx index of the source column.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  void multiply_source_and_add_to(const Field_element& coefficient, Index sourceColumnIndex, Index targetColumnIndex);

  /**
   * @brief Zeroes the entry at the given coordinates.
   *
   * @warning They will be no verification to ensure that the zeroing makes sense for the validity of a
   * boundary matrix of a filtered complex. So should be used while knowing what one is doing.
   *
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   */
  void zero_entry(Index columnIndex, Index rowIndex);
  /**
   * @brief Zeroes the column at the given index.
   *
   * @warning They will be no verification to ensure that the zeroing makes sense for the validity of a
   * boundary matrix of a filtered complex. So should be used while knowing what one is doing.
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
   * @brief Returns the pivot of the given column.
   *
   * @param columnIndex @ref MatIdx index of the column.
   * @return Pivot of the column at @p columnIndex.
   */
  Index get_pivot(Index columnIndex);

  /**
   * @brief Resets the matrix to an empty matrix.
   *
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  void reset(Column_settings* colSettings)
  {
    if constexpr (activeDimOption_) Dim_opt::_reset();
    if constexpr (activeSwapOption_) Swap_opt::_reset();
    if constexpr (activePairingOption_) Pair_opt::_reset();
    matrix_.clear();
    nextInsertIndex_ = 0;
    colSettings_ = colSettings;
  }

  /**
   * @brief Assign operator.
   */
  Boundary_matrix& operator=(const Boundary_matrix& other);
  /**
   * @brief Move assign operator.
   */
  Boundary_matrix& operator=(Boundary_matrix&& other) noexcept;

  /**
   * @brief Swap operator.
   */
  friend void swap(Boundary_matrix& matrix1, Boundary_matrix& matrix2) noexcept
  {
    swap(static_cast<Dim_opt&>(matrix1), static_cast<Dim_opt&>(matrix2));
    swap(static_cast<Swap_opt&>(matrix1), static_cast<Swap_opt&>(matrix2));
    swap(static_cast<Pair_opt&>(matrix1), static_cast<Pair_opt&>(matrix2));
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

  friend Swap_opt;
  friend Pair_opt;

  Column_container matrix_;      /**< Column container. */
  Index nextInsertIndex_;        /**< Next unused column index. */
  Column_settings* colSettings_; /**< Entry factory. */

  void _orderRowsIfNecessary();
  const Column& _get_column(Index columnIndex) const;
  Column& _get_column(Index columnIndex);
  Index _get_real_row_index(Index rowIndex) const;
  template <class Container>
  void _container_insert(const Container& column, Index pos, Dimension dim);
  void _container_insert(const Column& column, [[maybe_unused]] Index pos = 0);
};

template <class Master_matrix>
inline Boundary_matrix<Master_matrix>::Boundary_matrix(Column_settings* colSettings)
    : Dim_opt(Master_matrix::template get_null_value<Dimension>()),
      Swap_opt(),
      Pair_opt(),
      RA_opt(),
      nextInsertIndex_(0),
      colSettings_(colSettings)
{}

template <class Master_matrix>
template <class Boundary_range>
inline Boundary_matrix<Master_matrix>::Boundary_matrix(const std::vector<Boundary_range>& orderedBoundaries,
                                                       Column_settings* colSettings)
    : Dim_opt(Master_matrix::template get_null_value<Dimension>()),
      Swap_opt(orderedBoundaries.size()),
      Pair_opt(),
      RA_opt(orderedBoundaries.size()),
      nextInsertIndex_(orderedBoundaries.size()),
      colSettings_(colSettings)
{
  matrix_.reserve(orderedBoundaries.size());

  for (Index i = 0; i < orderedBoundaries.size(); i++) {
    _container_insert(orderedBoundaries[i], i, orderedBoundaries[i].size() == 0 ? 0 : orderedBoundaries[i].size() - 1);
  }
}

template <class Master_matrix>
inline Boundary_matrix<Master_matrix>::Boundary_matrix(unsigned int numberOfColumns, Column_settings* colSettings)
    : Dim_opt(Master_matrix::template get_null_value<Dimension>()),
      Swap_opt(numberOfColumns),
      Pair_opt(),
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
inline Boundary_matrix<Master_matrix>::Boundary_matrix(const Boundary_matrix& matrixToCopy,
                                                       Column_settings* colSettings)
    : Dim_opt(static_cast<const Dim_opt&>(matrixToCopy)),
      Swap_opt(static_cast<const Swap_opt&>(matrixToCopy)),
      Pair_opt(static_cast<const Pair_opt&>(matrixToCopy)),
      RA_opt(static_cast<const RA_opt&>(matrixToCopy)),
      nextInsertIndex_(matrixToCopy.nextInsertIndex_),
      colSettings_(colSettings == nullptr ? matrixToCopy.colSettings_ : colSettings)
{
  matrix_.reserve(matrixToCopy.matrix_.size());
  for (const auto& cont : matrixToCopy.matrix_) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      _container_insert(cont.second, cont.first);
    } else {
      _container_insert(cont);
    }
  }
}

template <class Master_matrix>
inline Boundary_matrix<Master_matrix>::Boundary_matrix(Boundary_matrix&& other) noexcept
    : Dim_opt(std::move(static_cast<Dim_opt&>(other))),
      Swap_opt(std::move(static_cast<Swap_opt&>(other))),
      Pair_opt(std::move(static_cast<Pair_opt&>(other))),
      RA_opt(std::move(static_cast<RA_opt&>(other))),
      matrix_(std::move(other.matrix_)),
      nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0)),
      colSettings_(std::exchange(other.colSettings_, nullptr))
{
}

template <class Master_matrix>
template <class Boundary_range>
inline typename Boundary_matrix<Master_matrix>::Index Boundary_matrix<Master_matrix>::insert_boundary(
    const Boundary_range& boundary,
    Dimension dim)
{
  return insert_boundary(nextInsertIndex_, boundary, dim);
}

template <class Master_matrix>
template <class Boundary_range>
inline typename Boundary_matrix<Master_matrix>::Index
Boundary_matrix<Master_matrix>::insert_boundary(ID_index cellIndex, const Boundary_range& boundary, Dimension dim)
{
  if (dim == Master_matrix::template get_null_value<Dimension>()) dim = boundary.size() == 0 ? 0 : boundary.size() - 1;

  _orderRowsIfNecessary();

  // updates container sizes
  if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_removable_rows) {
    if (boundary.size() != 0) {
      // row container
      RA_opt::_resize(Master_matrix::get_row_index(*std::prev(boundary.end())));
    }
  }

  // row swap map containers
  if constexpr (activeSwapOption_) {
    Swap_opt::_initialize_row_index(cellIndex);
  }

  // maps for possible shifting between column content and position indices used for birth events
  if constexpr (activePairingOption_) {
    if (cellIndex != nextInsertIndex_) {
      Pair_opt::_insert_id_position(cellIndex, nextInsertIndex_);
      if constexpr (Master_matrix::Option_list::has_removable_columns) {
        Pair_opt::PIDM::map_.emplace(nextInsertIndex_, cellIndex);
      }
    }
  }

  _container_insert(boundary, nextInsertIndex_, dim);

  return nextInsertIndex_++;
}

template <class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Column& Boundary_matrix<Master_matrix>::get_column(Index columnIndex)
{
  _orderRowsIfNecessary();

  return _get_column(columnIndex);
}

template <class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Row& Boundary_matrix<Master_matrix>::get_row(Index rowIndex)
{
  static_assert(Master_matrix::Option_list::has_row_access, "'get_row' is not implemented for the chosen options.");

  _orderRowsIfNecessary();

  return RA_opt::get_row(rowIndex);
}

template <class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Index Boundary_matrix<Master_matrix>::remove_last()
{
  static_assert(Master_matrix::Option_list::has_removable_columns,
                "'remove_last' is not implemented for the chosen options.");

  if (nextInsertIndex_ == 0) return Master_matrix::template get_null_value<Index>();  // empty matrix
  --nextInsertIndex_;

  // updates dimension max
  if constexpr (activeDimOption_) {
    Dim_opt::_update_down(matrix_.at(nextInsertIndex_).get_dimension());
  }

  // computes pivot and removes column from matrix_
  ID_index pivot;
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    auto it = matrix_.find(nextInsertIndex_);
    pivot = it->second.get_pivot();
    if constexpr (activeSwapOption_) {
      // if the removed column is positive, the pivot won't change value
      if (Swap_opt::_row_were_swapped() && pivot != Master_matrix::template get_null_value<ID_index>()) {
        Swap_opt::_orderRows();
        pivot = it->second.get_pivot();
      }
    }
    matrix_.erase(it);
  } else {
    pivot = matrix_[nextInsertIndex_].get_pivot();
    if constexpr (activeSwapOption_) {
      // if the removed column is positive, the pivot won't change value
      if (Swap_opt::_row_were_swapped() && pivot != Master_matrix::template get_null_value<ID_index>()) {
        Swap_opt::_orderRows();
        pivot = matrix_[nextInsertIndex_].get_pivot();
      }
    }
    if constexpr (Master_matrix::Option_list::has_row_access) {
      GUDHI_CHECK(nextInsertIndex_ == matrix_.size() - 1,
                  std::logic_error("Boundary_matrix::remove_last - Indexation problem."));
      matrix_.pop_back();
    } else {
      matrix_[nextInsertIndex_].clear();
    }
  }

  erase_empty_row(nextInsertIndex_);  // maximal, so empty

  // updates barcode
  if constexpr (activePairingOption_) {
    Pair_opt::_remove_last(nextInsertIndex_);
  }

  return pivot;
}

template <class Master_matrix>
inline void Boundary_matrix<Master_matrix>::erase_empty_row(Index rowIndex)
{
  // computes real row index and erases it if necessary from the row swap map containers
  ID_index rowID = rowIndex;
  if constexpr (activeSwapOption_) {
    rowID = Swap_opt::_erase_row(rowIndex);
  }

  if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows) {
    RA_opt::erase_empty_row(rowID);
  }
}

template <class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Index Boundary_matrix<Master_matrix>::get_number_of_columns() const
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.size();
  } else {
    return nextInsertIndex_;  // matrix could have been resized much bigger while insert
  }
}

template <class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Dimension Boundary_matrix<Master_matrix>::get_column_dimension(
    Index columnIndex) const
{
  return _get_column(columnIndex).get_dimension();
}

template <class Master_matrix>
inline void Boundary_matrix<Master_matrix>::add_to(Index sourceColumnIndex, Index targetColumnIndex)
{
  _get_column(targetColumnIndex) += _get_column(sourceColumnIndex);
}

template <class Master_matrix>
inline void Boundary_matrix<Master_matrix>::multiply_target_and_add_to(Index sourceColumnIndex,
                                                                       const Field_element& coefficient,
                                                                       Index targetColumnIndex)
{
  _get_column(targetColumnIndex).multiply_target_and_add(coefficient, _get_column(sourceColumnIndex));
}

template <class Master_matrix>
inline void Boundary_matrix<Master_matrix>::multiply_source_and_add_to(const Field_element& coefficient,
                                                                       Index sourceColumnIndex,
                                                                       Index targetColumnIndex)
{
  _get_column(targetColumnIndex).multiply_source_and_add(_get_column(sourceColumnIndex), coefficient);
}

template <class Master_matrix>
inline void Boundary_matrix<Master_matrix>::zero_entry(Index columnIndex, Index rowIndex)
{
  _get_column(columnIndex).clear(_get_real_row_index(rowIndex));
}

template <class Master_matrix>
inline void Boundary_matrix<Master_matrix>::zero_column(Index columnIndex)
{
  _get_column(columnIndex).clear();
}

template <class Master_matrix>
inline bool Boundary_matrix<Master_matrix>::is_zero_entry(Index columnIndex, Index rowIndex) const
{
  return !(_get_column(columnIndex).is_non_zero(_get_real_row_index(rowIndex)));
}

template <class Master_matrix>
inline bool Boundary_matrix<Master_matrix>::is_zero_column(Index columnIndex)
{
  return _get_column(columnIndex).is_empty();
}

template <class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Index Boundary_matrix<Master_matrix>::get_pivot(Index columnIndex)
{
  _orderRowsIfNecessary();

  return _get_column(columnIndex).get_pivot();
}

template <class Master_matrix>
inline Boundary_matrix<Master_matrix>& Boundary_matrix<Master_matrix>::operator=(const Boundary_matrix& other)
{
  if (this == &other) return *this;

  Dim_opt::operator=(other);
  Swap_opt::operator=(other);
  Pair_opt::operator=(other);
  RA_opt::operator=(other);

  matrix_.clear();
  nextInsertIndex_ = other.nextInsertIndex_;
  colSettings_ = other.colSettings_;

  matrix_.reserve(other.matrix_.size());
  for (const auto& cont : other.matrix_) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      _container_insert(cont.second, cont.first);
    } else {
      _container_insert(cont);
    }
  }

  return *this;
}

template <class Master_matrix>
inline Boundary_matrix<Master_matrix>& Boundary_matrix<Master_matrix>::operator=(Boundary_matrix&& other) noexcept
{
  if (this == &other) return *this;

  Dim_opt::operator=(std::move(other));
  Swap_opt::operator=(std::move(other));
  Pair_opt::operator=(std::move(other));
  RA_opt::operator=(std::move(other));

  matrix_ = std::move(other.matrix_);
  nextInsertIndex_ = std::exchange(other.nextInsertIndex_, 0);
  colSettings_ = std::exchange(other.colSettings_, nullptr);

  return *this;
}

template <class Master_matrix>
inline void Boundary_matrix<Master_matrix>::print()
{
  if constexpr (activeSwapOption_) {
    if (Swap_opt::_row_were_swapped()) Swap_opt::_orderRows();
  }
  std::cout << "Boundary_matrix:\n";
  for (Index i = 0; i < nextInsertIndex_; ++i) {
    Column& col = matrix_[i];
    for (auto e : col.get_content(nextInsertIndex_)) {
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
    for (ID_index i = 0; i < nextInsertIndex_; ++i) {
      const auto& row = RA_opt::get_row(i);
      for (const typename Column::Entry& entry : row) {
        std::cout << entry.get_column_index() << " ";
      }
      std::cout << "(" << i << ")\n";
    }
    std::cout << "\n";
  }
}

template <class Master_matrix>
inline void Boundary_matrix<Master_matrix>::_orderRowsIfNecessary()
{
  if constexpr (activeSwapOption_) {
    if (Swap_opt::_row_were_swapped()) Swap_opt::_orderRows();
  }
}

template <class Master_matrix>
inline const typename Boundary_matrix<Master_matrix>::Column& Boundary_matrix<Master_matrix>::_get_column(
    Index columnIndex) const
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.at(columnIndex);
  } else {
    return matrix_[columnIndex];
  }
}

template <class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Column& Boundary_matrix<Master_matrix>::_get_column(Index columnIndex)
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.at(columnIndex);
  } else {
    return matrix_[columnIndex];
  }
}

template <class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Index Boundary_matrix<Master_matrix>::_get_real_row_index(
    Index rowIndex) const
{
  if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update) {
    return Swap_opt::_get_row_index(rowIndex);
  } else {
    return rowIndex;
  }
}

template <class Master_matrix>
template <class Container>
inline void Boundary_matrix<Master_matrix>::_container_insert(const Container& column, Index pos, Dimension dim)
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
      if (matrix_.size() <= pos) {
        matrix_.emplace_back(column, dim, colSettings_);
      } else {
        matrix_[pos] = Column(column, dim, colSettings_);
      }
    }
  }
  if constexpr (activeDimOption_) {
    Dim_opt::_update_up(dim);
  }
}

template <class Master_matrix>
inline void Boundary_matrix<Master_matrix>::_container_insert(const Column& column, [[maybe_unused]] Index pos)
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.try_emplace(pos, Column(column, column.get_column_index(), RA_opt::_get_rows_ptr(), colSettings_));
    } else {
      matrix_.try_emplace(pos, Column(column, colSettings_));
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.emplace_back(column, column.get_column_index(), RA_opt::_get_rows_ptr(), colSettings_);
    } else {
      matrix_.emplace_back(column, colSettings_);
    }
  }
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_BOUNDARY_MATRIX_H
