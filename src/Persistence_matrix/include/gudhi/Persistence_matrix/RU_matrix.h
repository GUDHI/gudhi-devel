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
 * @file RU_matrix.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::RU_matrix class.
 */

#ifndef PM_RU_MATRIX_H
#define PM_RU_MATRIX_H

#include <utility>      //std::swap, std::move & std::exchange
#include <iostream>     //print() only
#include <vector>
#include <unordered_map>

namespace Gudhi {
namespace persistence_matrix {

template <class Master_matrix>
class RU_pairing;

/**
 * @class RU_matrix RU_matrix.h gudhi/Persistence_matrix/RU_matrix.h
 * @ingroup persistence_matrix
 *
 * @brief %Matrix structure to store the ordered @ref boundarymatrix "boundary matrix" \f$ R \cdot U \f$ of a filtered
 * complex in order to compute its persistent homology, as well as representative cycles.
 * Supports vineyards (see @cite vineyards) and the removal of maximal cells while maintaining
 * a valid barcode. Provides an access to its columns and rows.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class RU_matrix : public Master_matrix::RU_pairing_option,
                  public Master_matrix::RU_vine_swap_option,
                  public Master_matrix::RU_representative_cycles_option
{
 private:
  using Pair_opt = typename Master_matrix::RU_pairing_option;
  using Swap_opt = typename Master_matrix::RU_vine_swap_option;
  using Rep_opt = typename Master_matrix::RU_representative_cycles_option;

 public:
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element = typename Master_matrix::Element;               /**< Type of an field element. */
  using Column = typename Master_matrix::Column;                       /**< Column type. */
  using Row = typename Master_matrix::Row;                             /**< Row type,
                                                                            only necessary with row access option. */
  using Entry_constructor = typename Master_matrix::Entry_constructor; /**< Factory of @ref Entry classes. */
  using Column_settings = typename Master_matrix::Column_settings;     /**< Structure giving access to the columns to
                                                                            necessary external classes. */
  using Boundary = typename Master_matrix::Boundary;                   /**< Type of an input column. */
  using Index = typename Master_matrix::Index;                         /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;                   /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index;                 /**< @ref PosIdx index type. */
  using Dimension = typename Master_matrix::Dimension;                 /**< Dimension value type. */

  /**
   * @brief Constructs an empty matrix.
   *
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  RU_matrix(Column_settings* colSettings);
  /**
   * @brief Constructs a new matrix from the given ranges of @ref Matrix::Entry_representative. Each range corresponds
   * to a column (the order of the ranges are preserved). The content of the ranges is assumed to be sorted by
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
   * interest and not everything should be stored, then use the @ref insert_boundary method instead (after creating the
   * matrix with the @ref RU_matrix(unsigned int, Column_settings*) constructor preferably).
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  template <class Boundary_range = Boundary>
  RU_matrix(const std::vector<Boundary_range>& orderedBoundaries, Column_settings* colSettings);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   *
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the chosen column type, such as custom allocators.
   */
  RU_matrix(unsigned int numberOfColumns, Column_settings* colSettings);
  /**
   * @brief Copy constructor. If @p colSettings is not a null pointer, its value is kept
   * instead of the one in the copied matrix.
   *
   * @param matrixToCopy Matrix to copy.
   * @param colSettings Either a pointer to an existing setting structure for the columns or a null pointer.
   * The structure should contain all the necessary external classes specifically necessary for the chosen column type,
   * such as custom allocators. If null pointer, the pointer stored in @p matrixToCopy is used instead.
   */
  RU_matrix(const RU_matrix& matrixToCopy, Column_settings* colSettings = nullptr);
  /**
   * @brief Move constructor.
   *
   * @param other Matrix to move.
   */
  RU_matrix(RU_matrix&& other) noexcept;

  ~RU_matrix() = default;

  /**
   * @brief Inserts at the end of the matrix a new ordered column corresponding to the given boundary.
   * This means that it is assumed that this method is called on boundaries in the order of the filtration.
   * It also assumes that the cells in the given boundary are identified by their relative position in the filtration,
   * starting at 0. If it is not the case, use the other
   * @ref insert_boundary(ID_index, const Boundary_range&, Dimension) "insert_boundary" instead by indicating the
   * cell ID used in the boundaries when the cell is inserted.
   *
   * Different to the constructor, the boundaries do not have to come from a simplicial complex, but also from
   * a more general entry complex. This includes cubical complexes or Morse complexes for example.
   *
   * At the insertion, the boundary is stored in its reduced form and the barcode, if enabled, is updated.
   *
   * @tparam Boundary_range Range of @ref Matrix::Entry_representative. Assumed to have a begin(), end() and size()
   * method.
   * @param boundary Boundary generating the new column. The content should be ordered by ID.
   * @param dim Dimension of the cell whose boundary is given. If the complex is simplicial,
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   */
  template <class Boundary_range = Boundary>
  void insert_boundary(const Boundary_range& boundary,
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
   */
  template <class Boundary_range = Boundary>
  void insert_boundary(ID_index cellIndex,
                       const Boundary_range& boundary,
                       Dimension dim = Master_matrix::template get_null_value<Dimension>());
  /**
   * @brief Returns the column at the given @ref MatIdx index in \f$ R \f$ if @p inR is true and
   * in \f$ U \f$ if @p inR is false.
   * The type of the column depends on the chosen options, see @ref PersistenceMatrixOptions::column_type.
   *
   * Note that before returning the column, all column entries can eventually be reordered, if lazy swaps occurred.
   * It is therefore recommended to avoid calling @ref get_column between vine swaps, otherwise the benefits
   * of the the laziness is lost.
   *
   * @param columnIndex @ref MatIdx index of the column to return.
   * @param inR If true, returns the column in \f$ R \f$, if false, returns the column in \f$ U \f$.
   * Default value: true.
   * @return Reference to the column.
   */
  Column& get_column(Index columnIndex, bool inR = true);
  /**
   * @brief Returns the row at the given @ref rowindex "row index" in \f$ R \f$ if @p inR is true and
   * in \f$ U \f$ if @p inR is false.
   * The type of the row depends on the chosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   *
   * Note that before returning the row, all column entries can eventually be reordered, if lazy swaps occurred.
   * It is therefore recommended to avoid calling @ref get_row between vine swaps, otherwise the benefits
   * of the the laziness is lost.
   *
   * @param rowIndex @ref rowindex "Row index" of the row to return.
   * @param inR If true, returns the row in \f$ R \f$, if false, returns the row in \f$ U \f$.
   * Default value: true.
   * @return Reference to the row.
   */
  Row& get_row(Index rowIndex, bool inR = true);
  /**
   * @brief If @ref PersistenceMatrixOptions::has_row_access and @ref PersistenceMatrixOptions::has_removable_rows
   * are true: assumes that the row is empty in \f$ R \f$ and removes it from \f$ R \f$. If the matrix is valid,
   * a row will never be empty in \f$ U \f$, so \f$ U \f$ won't be affected.
   * If @ref PersistenceMatrixOptions::has_map_column_container
   * and @ref PersistenceMatrixOptions::has_column_and_row_swaps are true: cleans up maps used for the lazy row swaps.
   * Otherwise, does nothing.
   *
   * @warning The removed rows are always assumed to be empty in \f$ R \f$. If it is not the case, the deleted row
   * entries are not removed from their columns. And in the case of intrusive rows, this will generate a segmentation
   * fault when the column entries are destroyed later. The row access is just meant as a "read only" access to the
   * rows and the @ref erase_empty_row method just as a way to specify that a row is empty and can therefore be removed
   * from dictionaries. This allows to avoid testing the emptiness of a row at each column entry removal, what can
   * be quite frequent.
   *
   * @param rowIndex @ref rowindex "Row index" of the empty row.
   */
  void erase_empty_row(Index rowIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns and
   * @ref PersistenceMatrixOptions::has_vine_update are true.
   * Assumes that the cell is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., RU is still an upper triangular decomposition of the @ref boundarymatrix "boundary matrix").
   * The maximality of the cell is not verified.
   * Also updates the barcode if it is stored.
   *
   * See also @ref remove_last.
   *
   * @param columnIndex @ref MatIdx index of the cell to remove.
   */
  void remove_maximal_cell(Index columnIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns is true.
   * Removes the last cell in the filtration from the matrix and updates the barcode if it is stored.
   *
   * See also @ref remove_maximal_cell.
   */
  void remove_last();

  /**
   * @brief Returns the maximal dimension of a cell stored in the matrix.
   * Only available if @ref PersistenceMatrixOptions::has_matrix_maximal_dimension_access is true.
   *
   * @return The maximal dimension.
   */
  Dimension get_max_dimension() const;
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
   * @brief Zeroes the entry at the given coordinates in \f$ R \f$ if @p inR is true or in
   * \f$ U \f$ if @p inR is false. Should be used with care to not destroy the validity of the persistence
   * related properties of the matrix.
   *
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @param inR Boolean indicating in which matrix to zero: if true in \f$ R \f$ and if false in \f$ U \f$.
   * Default value: true.
   */
  void zero_entry(Index columnIndex, Index rowIndex, bool inR = true);
  /**
   * @brief Zeroes the column at the given index in \f$ R \f$ if @p inR is true or in
   * \f$ U \f$ if @p inR is false. Should be used with care to not destroy the validity of the persistence
   * related properties of the matrix.
   *
   * @param columnIndex @ref MatIdx index of the column to zero.
   * @param inR Boolean indicating in which matrix to zero: if true in \f$ R \f$ and if false in \f$ U \f$.
   * Default value: true.
   */
  void zero_column(Index columnIndex, bool inR = true);
  /**
   * @brief Indicates if the entry at given coordinates has value zero in \f$ R \f$
   * if @p inR is true or in \f$ U \f$ if @p inR is false.
   *
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @param inR Boolean indicating in which matrix to look: if true in \f$ R \f$ and if false in \f$ U \f$.
   * Default value: true.
   * @return true If the entry has value zero.
   * @return false Otherwise.
   */
  bool is_zero_entry(Index columnIndex, Index rowIndex, bool inR = true) const;
  /**
   * @brief Indicates if the column at given index has value zero in \f$ R \f$
   * if @p inR is true or in \f$ U \f$ if @p inR is false.
   *
   * Note that if @p inR is false, this method should usually return false.
   *
   * @param columnIndex @ref MatIdx index of the column.
   * @param inR Boolean indicating in which matrix to look: if true in \f$ R \f$ and if false in \f$ U \f$.
   * Default value: true.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(Index columnIndex, bool inR = true);

  /**
   * @brief Returns the @ref MatIdx index of the column which has the given @ref rowindex "row index" as pivot in
   * \f$ R \f$. Assumes that the pivot exists.
   *
   * @param cellIndex @ref rowindex "Row index" of the pivot.
   * @return @ref MatIdx index of the column in \f$ R \f$ with the given pivot.
   */
  Index get_column_with_pivot(Index cellIndex) const;
  /**
   * @brief Returns the @ref rowindex "row index" of the pivot of the given column in \f$ R \f$.
   *
   * @param columnIndex @ref MatIdx index of the column in \f$ R \f$.
   * @return The @ref rowindex "row index" of the pivot.
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
    if constexpr (Master_matrix::Option_list::has_column_pairings) Pair_opt::_reset();
    if constexpr (Master_matrix::Option_list::can_retrieve_representative_cycles) Rep_opt::_reset();
    reducedMatrixR_.reset(colSettings);
    mirrorMatrixU_.reset(colSettings);
    pivotToColumnIndex_.clear();
    nextEventIndex_ = 0;
    positionToID_.clear();
    operators_ = Master_matrix::get_operator_ptr(colSettings);
  }

  /**
   * @brief Assign operator.
   */
  RU_matrix& operator=(const RU_matrix& other);
  /**
   * @brief Move assign operator.
   */
  RU_matrix& operator=(RU_matrix&& other) noexcept;

  /**
   * @brief Swap operator.
   */
  friend void swap(RU_matrix& matrix1, RU_matrix& matrix2) noexcept
  {
    swap(static_cast<Pair_opt&>(matrix1), static_cast<Pair_opt&>(matrix2));
    swap(static_cast<Swap_opt&>(matrix1), static_cast<Swap_opt&>(matrix2));
    swap(static_cast<Rep_opt&>(matrix1), static_cast<Rep_opt&>(matrix2));
    swap(matrix1.reducedMatrixR_, matrix2.reducedMatrixR_);
    swap(matrix1.mirrorMatrixU_, matrix2.mirrorMatrixU_);
    matrix1.pivotToColumnIndex_.swap(matrix2.pivotToColumnIndex_);
    std::swap(matrix1.nextEventIndex_, matrix2.nextEventIndex_);
    matrix1.positionToID_.swap(matrix2.positionToID_);
    std::swap(matrix1.operators_, matrix2.operators_);
  }

  void print();  // for debug

 private:
  using Pivot_dictionary = typename Master_matrix::template Dictionary<Index>;
  using Position_dictionary = std::unordered_map<Pos_index, ID_index>;  // TODO: try other type of maps?
  using Barcode = typename Master_matrix::Barcode;
  using Bar_dictionary = typename Master_matrix::Bar_dictionary;
  using R_matrix = typename Master_matrix::Master_boundary_matrix;
  using U_matrix = typename Master_matrix::Master_base_matrix;

  friend Rep_opt;                    // direct access to the two matrices
  friend Swap_opt;                   // direct access to the two matrices, pivotToColumnIndex_
  friend RU_pairing<Master_matrix>;  // direct access to positionToID_

  R_matrix reducedMatrixR_; /**< R. */
  // TODO: make U not accessible by default and add option to enable access? Inaccessible, it
  // needs less options and we could avoid some ifs.
  U_matrix mirrorMatrixU_;              /**< U. */
  Pivot_dictionary pivotToColumnIndex_; /**< Map from pivot row index to column @ref MatIdx index. */
  Pos_index nextEventIndex_;            /**< Next birth or death index. */
  Position_dictionary positionToID_;    /**< Map from @ref MatIdx to @ref IDIdx. */
  Field_operators const* operators_;    /**< Field operators, can be nullptr if
                                             @ref PersistenceMatrixOptions::is_z2 is true. */

  void _insert_boundary(Index currentIndex);
  void _initialize_U();
  void _reduce();
  void _reduce_last_column(Index lastIndex);
  void _reduce_column(Index target, Index eventIndex);
  void _reduce_column_by(Index target, Index source);
  Index _get_column_with_pivot(ID_index pivot) const;
  void _update_barcode(ID_index birthPivot, Pos_index death);
  void _add_bar(Dimension dim, Pos_index birth);
  void _remove_last_in_barcode(Pos_index eventIndex);
};

template <class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(Column_settings* colSettings)
    : Pair_opt(),
      Swap_opt(),
      Rep_opt(),
      reducedMatrixR_(colSettings),
      mirrorMatrixU_(colSettings),
      nextEventIndex_(0),
      operators_(Master_matrix::get_operator_ptr(colSettings))
{}

template <class Master_matrix>
template <class Boundary_range>
inline RU_matrix<Master_matrix>::RU_matrix(const std::vector<Boundary_range>& orderedBoundaries,
                                           Column_settings* colSettings)
    : Pair_opt(),
      Swap_opt(),
      Rep_opt(),
      reducedMatrixR_(orderedBoundaries, colSettings),
      mirrorMatrixU_(orderedBoundaries.size(), colSettings),
      nextEventIndex_(orderedBoundaries.size()),
      operators_(Master_matrix::get_operator_ptr(colSettings))
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.reserve(orderedBoundaries.size());
  } else {
    pivotToColumnIndex_.resize(orderedBoundaries.size(), Master_matrix::template get_null_value<Index>());
  }

  _initialize_U();
  _reduce();
}

template <class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(unsigned int numberOfColumns, Column_settings* colSettings)
    : Pair_opt(),
      Swap_opt(),
      Rep_opt(),
      reducedMatrixR_(numberOfColumns, colSettings),
      mirrorMatrixU_(numberOfColumns, colSettings),
      nextEventIndex_(0),
      positionToID_(numberOfColumns),
      operators_(Master_matrix::get_operator_ptr(colSettings))
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.reserve(numberOfColumns);
  } else {
    pivotToColumnIndex_.resize(numberOfColumns, Master_matrix::template get_null_value<Index>());
  }
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    Pair_opt::_reserve(numberOfColumns);
  }
}

template <class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(const RU_matrix& matrixToCopy, Column_settings* colSettings)
    : Pair_opt(static_cast<const Pair_opt&>(matrixToCopy)),
      Swap_opt(static_cast<const Swap_opt&>(matrixToCopy)),
      Rep_opt(static_cast<const Rep_opt&>(matrixToCopy)),
      reducedMatrixR_(matrixToCopy.reducedMatrixR_, colSettings),
      mirrorMatrixU_(matrixToCopy.mirrorMatrixU_, colSettings),
      pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
      nextEventIndex_(matrixToCopy.nextEventIndex_),
      positionToID_(matrixToCopy.positionToID_),
      operators_(colSettings == nullptr ? matrixToCopy.operators_ : Master_matrix::get_operator_ptr(colSettings))
{}

template <class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(RU_matrix&& other) noexcept
    : Pair_opt(std::move(static_cast<Pair_opt&>(other))),
      Swap_opt(std::move(static_cast<Swap_opt&>(other))),
      Rep_opt(std::move(static_cast<Rep_opt&>(other))),
      reducedMatrixR_(std::move(other.reducedMatrixR_)),
      mirrorMatrixU_(std::move(other.mirrorMatrixU_)),
      pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
      nextEventIndex_(std::exchange(other.nextEventIndex_, 0)),
      positionToID_(std::move(other.positionToID_)),
      operators_(std::exchange(other.operators_, nullptr))
{}

template <class Master_matrix>
template <class Boundary_range>
inline void RU_matrix<Master_matrix>::insert_boundary(const Boundary_range& boundary, Dimension dim)
{
  _insert_boundary(reducedMatrixR_.insert_boundary(boundary, dim));
}

template <class Master_matrix>
template <class Boundary_range>
inline void RU_matrix<Master_matrix>::insert_boundary(ID_index cellIndex, const Boundary_range& boundary, Dimension dim)
{
  // maps for possible shifting between column content and position indices used for birth events
  if (cellIndex != nextEventIndex_) {
    positionToID_.emplace(nextEventIndex_, cellIndex);
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      Pair_opt::_insert_id_position(cellIndex, nextEventIndex_);
    }
  }

  _insert_boundary(reducedMatrixR_.insert_boundary(cellIndex, boundary, dim));
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Column& RU_matrix<Master_matrix>::get_column(Index columnIndex, bool inR)
{
  if (inR) {
    return reducedMatrixR_.get_column(columnIndex);
  }
  return mirrorMatrixU_.get_column(columnIndex);
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Row& RU_matrix<Master_matrix>::get_row(Index rowIndex, bool inR)
{
  static_assert(Master_matrix::Option_list::has_row_access, "'get_row' is not implemented for the chosen options.");

  if (inR) {
    return reducedMatrixR_.get_row(rowIndex);
  }
  return mirrorMatrixU_.get_row(rowIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::erase_empty_row(Index rowIndex)
{
  reducedMatrixR_.erase_empty_row(rowIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::remove_maximal_cell(Index columnIndex)
{
  static_assert(Master_matrix::Option_list::has_removable_columns && Master_matrix::Option_list::has_vine_update,
                "'remove_maximal_cell' is not implemented for the chosen options.");

  // TODO: is there an easy test to verify maximality even without row access?

  for (Index curr = columnIndex; curr < nextEventIndex_ - 1; ++curr) {
    Swap_opt::vine_swap(curr);
  }

  remove_last();
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::remove_last()
{
  static_assert(Master_matrix::Option_list::has_removable_columns,
                "'remove_last' is not implemented for the chosen options.");

  if (nextEventIndex_ == 0) return;  // empty matrix
  --nextEventIndex_;

  // assumes PosIdx == MatIdx for boundary matrices.
  _remove_last_in_barcode(nextEventIndex_);

  mirrorMatrixU_.remove_last();
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.erase(reducedMatrixR_.remove_last());
  } else {
    ID_index lastPivot = reducedMatrixR_.remove_last();
    if (lastPivot != Master_matrix::template get_null_value<ID_index>())
      pivotToColumnIndex_[lastPivot] = Master_matrix::template get_null_value<Index>();
  }

  // if has_column_pairings is true, then the element is already removed in _remove_last_in_barcode
  // to avoid a second "find"
  if constexpr (!Master_matrix::Option_list::has_column_pairings) {
    positionToID_.erase(nextEventIndex_);
  }
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Dimension RU_matrix<Master_matrix>::get_max_dimension() const
{
  return reducedMatrixR_.get_max_dimension();
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Index RU_matrix<Master_matrix>::get_number_of_columns() const
{
  return reducedMatrixR_.get_number_of_columns();
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Dimension RU_matrix<Master_matrix>::get_column_dimension(
    Index columnIndex) const
{
  return reducedMatrixR_.get_column_dimension(columnIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(Index sourceColumnIndex, Index targetColumnIndex)
{
  reducedMatrixR_.add_to(sourceColumnIndex, targetColumnIndex);
  // U transposed to avoid row operations
  if constexpr (Master_matrix::Option_list::is_z2)
    mirrorMatrixU_.add_to(targetColumnIndex, sourceColumnIndex);
  else
    mirrorMatrixU_.multiply_source_and_add_to(
        operators_->get_characteristic() - 1, targetColumnIndex, sourceColumnIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::multiply_target_and_add_to(Index sourceColumnIndex,
                                                                 const Field_element& coefficient,
                                                                 Index targetColumnIndex)
{
  reducedMatrixR_.multiply_target_and_add_to(sourceColumnIndex, coefficient, targetColumnIndex);
  // U transposed to avoid row operations
  mirrorMatrixU_.get_column(targetColumnIndex) *= coefficient;
  mirrorMatrixU_.multiply_source_and_add_to(operators_->get_characteristic() - 1, targetColumnIndex, sourceColumnIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::multiply_source_and_add_to(const Field_element& coefficient,
                                                                 Index sourceColumnIndex,
                                                                 Index targetColumnIndex)
{
  reducedMatrixR_.multiply_source_and_add_to(coefficient, sourceColumnIndex, targetColumnIndex);
  // U transposed to avoid row operations
  if constexpr (Master_matrix::Option_list::is_z2) {
    if (coefficient) mirrorMatrixU_.add_to(targetColumnIndex, sourceColumnIndex);
  } else {
    mirrorMatrixU_.multiply_source_and_add_to(
        operators_->get_characteristic() - coefficient, targetColumnIndex, sourceColumnIndex);
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::zero_entry(Index columnIndex, Index rowIndex, bool inR)
{
  if (inR) {
    return reducedMatrixR_.zero_entry(columnIndex, rowIndex);
  }
  return mirrorMatrixU_.zero_entry(columnIndex, rowIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::zero_column(Index columnIndex, bool inR)
{
  if (inR) {
    return reducedMatrixR_.zero_column(columnIndex);
  }
  return mirrorMatrixU_.zero_column(columnIndex);
}

template <class Master_matrix>
inline bool RU_matrix<Master_matrix>::is_zero_entry(Index columnIndex, Index rowIndex, bool inR) const
{
  if (inR) {
    return reducedMatrixR_.is_zero_entry(columnIndex, rowIndex);
  }
  return mirrorMatrixU_.is_zero_entry(columnIndex, rowIndex);
}

template <class Master_matrix>
inline bool RU_matrix<Master_matrix>::is_zero_column(Index columnIndex, bool inR)
{
  if (inR) {
    return reducedMatrixR_.is_zero_column(columnIndex);
  }
  return mirrorMatrixU_.is_zero_column(columnIndex);
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Index RU_matrix<Master_matrix>::get_column_with_pivot(Index cellIndex) const
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return pivotToColumnIndex_.at(cellIndex);
  } else {
    return pivotToColumnIndex_[cellIndex];
  }
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Index RU_matrix<Master_matrix>::get_pivot(Index columnIndex)
{
  return reducedMatrixR_.get_column(columnIndex).get_pivot();
}

template <class Master_matrix>
inline RU_matrix<Master_matrix>& RU_matrix<Master_matrix>::operator=(const RU_matrix& other)
{
  if (this == &other) return *this;

  Swap_opt::operator=(other);
  Pair_opt::operator=(other);
  Rep_opt::operator=(other);
  reducedMatrixR_ = other.reducedMatrixR_;
  mirrorMatrixU_ = other.mirrorMatrixU_;
  pivotToColumnIndex_ = other.pivotToColumnIndex_;
  nextEventIndex_ = other.nextEventIndex_;
  positionToID_ = other.positionToID_;
  operators_ = other.operators_;

  return *this;
}

template <class Master_matrix>
inline RU_matrix<Master_matrix>& RU_matrix<Master_matrix>::operator=(RU_matrix&& other) noexcept
{
  if (this == &other) return *this;

  Pair_opt::operator=(std::move(other));
  Swap_opt::operator=(std::move(other));
  Rep_opt::operator=(std::move(other));

  reducedMatrixR_ = std::move(other.reducedMatrixR_);
  mirrorMatrixU_ = std::move(other.mirrorMatrixU_);
  pivotToColumnIndex_ = std::move(other.pivotToColumnIndex_);
  nextEventIndex_ = std::exchange(other.nextEventIndex_, 0);
  positionToID_ = std::move(other.positionToID_);
  operators_ = std::exchange(other.operators_, nullptr);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::print()
{
  std::cout << "R_matrix:\n";
  reducedMatrixR_.print();
  std::cout << "U_matrix:\n";
  mirrorMatrixU_.print();
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_insert_boundary(Index currentIndex)
{
  mirrorMatrixU_.insert_column(currentIndex, 1);

  if constexpr (!Master_matrix::Option_list::has_map_column_container) {
    ID_index pivot = reducedMatrixR_.get_column(currentIndex).get_pivot();
    if (pivot != Master_matrix::template get_null_value<ID_index>() && pivotToColumnIndex_.size() <= pivot)
      pivotToColumnIndex_.resize((pivot + 1) * 2, Master_matrix::template get_null_value<Index>());
  }

  _reduce_last_column(currentIndex);
  ++nextEventIndex_;
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_initialize_U()
{
  for (ID_index i = 0; i < reducedMatrixR_.get_number_of_columns(); i++) {
    mirrorMatrixU_.insert_column(i, 1);
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce()
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    Pair_opt::_reserve(reducedMatrixR_.get_number_of_columns());
  }

  for (Index i = 0; i < reducedMatrixR_.get_number_of_columns(); i++) {
    if (!(reducedMatrixR_.is_zero_column(i))) {
      _reduce_column(i, i);
    } else {
      _add_bar(get_column_dimension(i), i);
    }
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce_last_column(Index lastIndex)
{
  if (reducedMatrixR_.get_column(lastIndex).is_empty()) {
    _add_bar(get_column_dimension(lastIndex), nextEventIndex_);
    return;
  }

  _reduce_column(lastIndex, nextEventIndex_);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce_column(Index target, Index eventIndex)
{
  Column& curr = reducedMatrixR_.get_column(target);
  ID_index pivot = curr.get_pivot();
  Index currIndex = _get_column_with_pivot(pivot);

  while (pivot != Master_matrix::template get_null_value<ID_index>() &&
         currIndex != Master_matrix::template get_null_value<Index>()) {
    _reduce_column_by(target, currIndex);
    pivot = curr.get_pivot();
    currIndex = _get_column_with_pivot(pivot);
  }

  if (pivot != Master_matrix::template get_null_value<ID_index>()) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      pivotToColumnIndex_.try_emplace(pivot, target);
    } else {
      pivotToColumnIndex_[pivot] = target;
    }
    _update_barcode(pivot, eventIndex);
  } else {
    _add_bar(get_column_dimension(target), eventIndex);
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce_column_by(Index target, Index source)
{
  Column& curr = reducedMatrixR_.get_column(target);
  if constexpr (Master_matrix::Option_list::is_z2) {
    curr += reducedMatrixR_.get_column(source);
    // to avoid having to do line operations, U is transposed
    // TODO: explain this somewhere in the documentation...
    mirrorMatrixU_.get_column(source).push_back(*mirrorMatrixU_.get_column(target).begin());
  } else {
    Column& toadd = reducedMatrixR_.get_column(source);
    Field_element coef = toadd.get_pivot_value();
    coef = operators_->get_inverse(coef);
    operators_->multiply_inplace(coef, operators_->get_characteristic() - curr.get_pivot_value());

    curr.multiply_source_and_add(toadd, coef);
    auto entry = *mirrorMatrixU_.get_column(target).begin();
    operators_->multiply_inplace(entry.get_element(), operators_->get_characteristic() - coef);
    // to avoid having to do line operations, U is transposed
    // TODO: explain this somewhere in the documentation...
    mirrorMatrixU_.get_column(source).push_back(entry);
  }
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Index RU_matrix<Master_matrix>::_get_column_with_pivot(ID_index pivot) const
{
  if (pivot == Master_matrix::template get_null_value<ID_index>())
    return Master_matrix::template get_null_value<Index>();
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    auto it = pivotToColumnIndex_.find(pivot);
    if (it == pivotToColumnIndex_.end()) return Master_matrix::template get_null_value<Index>();
    return it->second;
  } else {
    if (pivot >= pivotToColumnIndex_.size()) return Master_matrix::template get_null_value<Index>();
    return pivotToColumnIndex_[pivot];
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_update_barcode(ID_index birthPivot, Pos_index death)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    Pair_opt::_update_barcode(birthPivot, death);
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_add_bar(Dimension dim, Pos_index birth)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    Pair_opt::_add_bar(dim, birth);
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_remove_last_in_barcode(Pos_index eventIndex)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    Pair_opt::_remove_last(eventIndex);
  }
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_RU_MATRIX_H
