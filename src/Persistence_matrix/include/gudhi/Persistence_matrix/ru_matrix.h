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
 * @file ru_matrix.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::RU_matrix class.
 */

#ifndef PM_RU_MATRIX_H
#define PM_RU_MATRIX_H

#include <vector>
#include <utility>   //std::swap, std::move & std::exchange
#include <iostream>  //print() only

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class RU_matrix ru_matrix.h gudhi/Persistence_matrix/ru_matrix.h
 * @ingroup persistence_matrix
 *
 * @brief %Matrix structure to store the ordered @ref boundarymatrix "boundary matrix" \f$ R \cdot U \f$ of a filtered
 * complex in order to compute its persistent homology, as well as representative cycles.
 * Supports vineyards (see @cite vineyards) and the removal of maximal faces while maintaining
 * a valid barcode. Provides an access to its columns and rows.
 * 
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class RU_matrix : public Master_matrix::RU_pairing_option,
                  public Master_matrix::RU_vine_swap_option,
                  public Master_matrix::RU_representative_cycles_option 
{
 public:
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element_type = typename Master_matrix::element_type;    /**< Type of an field element. */
  using Column_type = typename Master_matrix::Column_type;            /**< Column type. */
  using Row_type = typename Master_matrix::Row_type;                  /**< Row type,
                                                                           only necessary with row access option. */
  using Cell_constructor = typename Master_matrix::Cell_constructor;  /**< Factory of @ref Cell classes. */
  using Column_settings = typename Master_matrix::Column_settings;    /**< Structure giving access to the columns to
                                                                           necessary external classes. */
  using boundary_type = typename Master_matrix::boundary_type;        /**< Type of an input column. */
  using index = typename Master_matrix::index;                        /**< @ref MatIdx index type. */
  using id_index = typename Master_matrix::id_index;                  /**< @ref IDIdx index type. */
  using pos_index = typename Master_matrix::pos_index;                /**< @ref PosIdx index type. */
  using dimension_type = typename Master_matrix::dimension_type;      /**< Dimension value type. */

  /**
   * @brief Constructs an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  RU_matrix(Column_settings* colSettings);
  /**
   * @brief Constructs a new matrix from the given ranges of @ref Matrix::cell_rep_type. Each range corresponds to a
   * column (the order of the ranges are preserved). The content of the ranges is assumed to be sorted by increasing
   * IDs. The IDs of the simplices are also assumed to be consecutive, ordered by filtration value, starting with 0. 
   * 
   * @tparam Boundary_type Range type for @ref Matrix::cell_rep_type ranges.
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
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  template <class Boundary_type = boundary_type>
  RU_matrix(const std::vector<Boundary_type>& orderedBoundaries, 
            Column_settings* colSettings);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  RU_matrix(unsigned int numberOfColumns, Column_settings* colSettings);
  /**
   * @brief Copy constructor. If @p colSettings is not a null pointer, its value is kept
   * instead of the one in the copied matrix.
   * 
   * @param matrixToCopy Matrix to copy.
   * @param colSettings Either a pointer to an existing setting structure for the columns or a null pointer.
   * The structure should contain all the necessary external classes specifically necessary for the choosen column type,
   * such as custom allocators. If null pointer, the pointer stored in @p matrixToCopy is used instead.
   */
  RU_matrix(const RU_matrix& matrixToCopy, 
            Column_settings* colSettings = nullptr);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  RU_matrix(RU_matrix&& other) noexcept;

  /**
   * @brief Inserts at the end of the matrix a new ordered column corresponding to the given boundary.
   * This means that it is assumed that this method is called on boundaries in the order of the filtration.
   * It also assumes that the faces in the given boundary are identified by their relative position in the filtration,
   * starting at 0. If it is not the case, use the other
   * @ref insert_boundary(id_index, const Boundary_type&, dimension_type) "insert_boundary" instead by indicating the
   * face ID used in the boundaries when the face is inserted.
   *
   * Different to the constructor, the boundaries do not have to come from a simplicial complex, but also from
   * a more general cell complex. This includes cubical complexes or Morse complexes for example.
   *
   * At the insertion, the boundary is stored in its reduced form and the barcode, if enabled, is updated.
   * 
   * @tparam Boundary_type Range of @ref Matrix::cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param boundary Boundary generating the new column. The content should be ordered by ID.
   * @param dim Dimension of the face whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   */
  template <class Boundary_type = boundary_type>
  void insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
  /**
   * @brief It does the same as the other version, but allows the boundary faces to be identified without restrictions
   * except that all IDs have to be strictly increasing in the order of filtration. Note that you should avoid then
   * to use the other insertion method to avoid overwriting IDs.
   *
   * As a face has to be inserted before one of its cofaces in a valid filtration (recall that it is assumed that
   * the faces are inserted by order of filtration), it is sufficient to indicate the ID of the face being inserted.
   * 
   * @tparam Boundary_type Range of @ref Matrix::cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param faceIndex @ref IDIdx index to use to identify the new face.
   * @param boundary Boundary generating the new column. The indices of the boundary have to correspond to the 
   * @p faceIndex values of precedent calls of the method for the corresponding faces and should be ordered in 
   * increasing order.
   * @param dim Dimension of the face whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   */
  template <class Boundary_type = boundary_type>
  void insert_boundary(id_index faceIndex, const Boundary_type& boundary, dimension_type dim = -1);
  /**
   * @brief Returns the column at the given @ref MatIdx index in \f$ R \f$ if @p inR is true and
   * in \f$ U \f$ if @p inR is false.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   *
   * Note that before returning the column, all column cells can eventually be reordered, if lazy swaps occurred.
   * It is therefore recommended to avoid calling @ref get_column between vine swaps, otherwise the benefits
   * of the the laziness is lost.
   * 
   * @param columnIndex @ref MatIdx index of the column to return.
   * @param inR If true, returns the column in \f$ R \f$, if false, returns the column in \f$ U \f$.
   * Default value: true.
   * @return Reference to the column.
   */
  Column_type& get_column(index columnIndex, bool inR = true);
  /**
   * @brief Returns the row at the given @ref rowindex "row index" in \f$ R \f$ if @p inR is true and
   * in \f$ U \f$ if @p inR is false.
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   *
   * Note that before returning the row, all column cells can eventually be reordered, if lazy swaps occurred.
   * It is therefore recommended to avoid calling @ref get_row between vine swaps, otherwise the benefits
   * of the the laziness is lost.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return.
   * @param inR If true, returns the row in \f$ R \f$, if false, returns the row in \f$ U \f$.
   * Default value: true.
   * @return Reference to the row.
   */
  Row_type& get_row(index rowIndex, bool inR = true);
  /**
   * @brief If @ref PersistenceMatrixOptions::has_row_access and @ref PersistenceMatrixOptions::has_removable_rows
   * are true: assumes that the row is empty in \f$ R \f$ and removes it from \f$ R \f$. If the matrix is valid,
   * a row will never be empty in \f$ U \f$, so \f$ U \f$ won't be affected.
   * If @ref PersistenceMatrixOptions::has_map_column_container
   * and @ref PersistenceMatrixOptions::has_column_and_row_swaps are true: cleans up maps used for the lazy row swaps.
   * Otherwise, does nothing.
   *
   * @warning The removed rows are always assumed to be empty in \f$ R \f$. If it is not the case, the deleted row
   * cells are not removed from their columns. And in the case of intrusive rows, this will generate a segmentation
   * fault when the column cells are destroyed later. The row access is just meant as a "read only" access to the
   * rows and the @ref erase_empty_row method just as a way to specify that a row is empty and can therefore be removed
   * from dictionaries. This allows to avoid testing the emptiness of a row at each column cell removal, what can
   * be quite frequent. 
   * 
   * @param rowIndex @ref rowindex "Row index" of the empty row.
   */
  void erase_empty_row(index rowIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns and
   * @ref PersistenceMatrixOptions::has_vine_update are true.
   * Assumes that the face is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., RU is still an upper triangular decomposition of the @ref boundarymatrix "boundary matrix").
   * The maximality of the face is not verified.
   * Also updates the barcode if it is stored.
   *
   * See also @ref remove_last.
   * 
   * @param columnIndex @ref MatIdx index of the face to remove.
   */
  void remove_maximal_face(index columnIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns is true.
   * Removes the last face in the filtration from the matrix and updates the barcode if it is stored.
   *
   * See also @ref remove_maximal_face.
   */
  void remove_last();

  /**
   * @brief Returns the maximal dimension of a face stored in the matrix.
   * Only available if @ref PersistenceMatrixOptions::has_matrix_maximal_dimension_access is true.
   * 
   * @return The maximal dimension.
   */
  dimension_type get_max_dimension() const;
  /**
   * @brief Returns the current number of columns in the matrix.
   * 
   * @return The number of columns.
   */
  index get_number_of_columns() const;
  /**
   * @brief Returns the dimension of the given column.
   * 
   * @param columnIndex @ref MatIdx index of the column representing the face.
   * @return Dimension of the face.
   */
  dimension_type get_column_dimension(index columnIndex) const;

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
  void add_to(index sourceColumnIndex, index targetColumnIndex);
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
  void multiply_target_and_add_to(index sourceColumnIndex, 
                                  const Field_element_type& coefficient,
                                  index targetColumnIndex);
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
  void multiply_source_and_add_to(const Field_element_type& coefficient, 
                                  index sourceColumnIndex,
                                  index targetColumnIndex);

  /**
   * @brief Zeroes the cell at the given coordinates in \f$ R \f$ if @p inR is true or in
   * \f$ U \f$ if @p inR is false. Should be used with care to not destroy the validity of the persistence
   * related properties of the matrix.
   * 
   * @param columnIndex @ref MatIdx index of the column of the cell.
   * @param rowIndex @ref rowindex "Row index" of the row of the cell.
   * @param inR Boolean indicating in which matrix to zero: if true in \f$ R \f$ and if false in \f$ U \f$.
   * Default value: true.
   */
  void zero_cell(index columnIndex, index rowIndex, bool inR = true);
  /**
   * @brief Zeroes the column at the given index in \f$ R \f$ if @p inR is true or in
   * \f$ U \f$ if @p inR is false. Should be used with care to not destroy the validity of the persistence
   * related properties of the matrix.
   * 
   * @param columnIndex @ref MatIdx index of the column to zero.
   * @param inR Boolean indicating in which matrix to zero: if true in \f$ R \f$ and if false in \f$ U \f$.
   * Default value: true.
   */
  void zero_column(index columnIndex, bool inR = true);
  /**
   * @brief Indicates if the cell at given coordinates has value zero in \f$ R \f$
   * if @p inR is true or in \f$ U \f$ if @p inR is false.
   * 
   * @param columnIndex @ref MatIdx index of the column of the cell.
   * @param rowIndex @ref rowindex "Row index" of the row of the cell.
   * @param inR Boolean indicating in which matrix to look: if true in \f$ R \f$ and if false in \f$ U \f$.
   * Default value: true.
   * @return true If the cell has value zero.
   * @return false Otherwise.
   */
  bool is_zero_cell(index columnIndex, index rowIndex, bool inR = true) const;
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
  bool is_zero_column(index columnIndex, bool inR = true);

  /**
   * @brief Returns the @ref MatIdx index of the column which has the given @ref rowindex "row index" as pivot in
   * \f$ R \f$. Assumes that the pivot exists.
   * 
   * @param faceIndex @ref rowindex "Row index" of the pivot.
   * @return @ref MatIdx index of the column in \f$ R \f$ with the given pivot.
   */
  index get_column_with_pivot(index faceIndex) const;
  /**
   * @brief Returns the @ref rowindex "row index" of the pivot of the given column in \f$ R \f$.
   * 
   * @param columnIndex @ref MatIdx index of the column in \f$ R \f$.
   * @return The @ref rowindex "row index" of the pivot.
   */
  index get_pivot(index columnIndex);

  /**
   * @brief Resets the matrix to an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  void reset(Column_settings* colSettings) {
    reducedMatrixR_.reset(colSettings);
    mirrorMatrixU_.reset(colSettings);
    pivotToColumnIndex_.clear();
    nextEventIndex_ = 0;
    if constexpr (!Master_matrix::Option_list::is_z2){
      operators_ = &(colSettings->operators);
    }
  }

  /**
   * @brief Assign operator.
   */
  RU_matrix& operator=(const RU_matrix& other);
  /**
   * @brief Swap operator.
   */
  friend void swap(RU_matrix& matrix1, RU_matrix& matrix2) {
    swap(static_cast<typename Master_matrix::RU_pairing_option&>(matrix1),
         static_cast<typename Master_matrix::RU_pairing_option&>(matrix2));
    swap(static_cast<typename Master_matrix::RU_vine_swap_option&>(matrix1),
         static_cast<typename Master_matrix::RU_vine_swap_option&>(matrix2));
    swap(static_cast<typename Master_matrix::RU_representative_cycles_option&>(matrix1),
         static_cast<typename Master_matrix::RU_representative_cycles_option&>(matrix2));
    swap(matrix1.reducedMatrixR_, matrix2.reducedMatrixR_);
    swap(matrix1.mirrorMatrixU_, matrix2.mirrorMatrixU_);
    matrix1.pivotToColumnIndex_.swap(matrix2.pivotToColumnIndex_);
    std::swap(matrix1.nextEventIndex_, matrix2.nextEventIndex_);
    std::swap(matrix1.operators_, matrix2.operators_);
  }

  void print();  // for debug

 private:
  using swap_opt = typename Master_matrix::RU_vine_swap_option;
  using pair_opt = typename Master_matrix::RU_pairing_option;
  using rep_opt = typename Master_matrix::RU_representative_cycles_option;
  using dictionary_type = typename Master_matrix::template dictionary_type<index>;
  using barcode_type = typename Master_matrix::barcode_type;
  using bar_dictionary_type = typename Master_matrix::bar_dictionary_type;
  using r_matrix_type = typename Master_matrix::Boundary_matrix_type;
  using u_matrix_type = typename Master_matrix::Base_matrix_type;

  friend rep_opt;   // direct access to the two matrices
  friend swap_opt;  // direct access to the two matrices

  r_matrix_type reducedMatrixR_;        /**< R. */
  // TODO: make U not accessible by default and add option to enable access? Inaccessible, it
  // needs less options and we could avoid some ifs.
  u_matrix_type mirrorMatrixU_;         /**< U. */
  dictionary_type pivotToColumnIndex_; /**< Map from pivot row index to column @ref MatIdx index. */
  pos_index nextEventIndex_;            /**< Next birth or death index. */
  Field_operators* operators_;          /**< Field operators,
                                             can be nullptr if @ref PersistenceMatrixOptions::is_z2 is true. */

  void _insert_boundary(index currentIndex);
  void _initialize_U();
  void _reduce();
  void _reduce_last_column(index lastIndex);
  void _reduce_column(index target, index eventIndex);
  void _reduce_column_by(index target, index source);
  void _update_barcode(id_index birthPivot, pos_index death);
  void _add_bar(dimension_type dim, pos_index birth);
  void _remove_last_in_barcode(pos_index eventIndex);

  constexpr bar_dictionary_type& _indexToBar();
};

template <class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(Column_settings* colSettings)
    : pair_opt(),
      swap_opt(),
      rep_opt(),
      reducedMatrixR_(colSettings),
      mirrorMatrixU_(colSettings),
      nextEventIndex_(0),
      operators_(nullptr) 
{
  if constexpr (!Master_matrix::Option_list::is_z2){
    operators_ = &(colSettings->operators);
  }
}

template <class Master_matrix>
template <class Boundary_type>
inline RU_matrix<Master_matrix>::RU_matrix(const std::vector<Boundary_type>& orderedBoundaries,
                                           Column_settings* colSettings)
    : pair_opt(),
      swap_opt(),
      rep_opt(),
      reducedMatrixR_(orderedBoundaries, colSettings),
      mirrorMatrixU_(orderedBoundaries.size(), colSettings),
      nextEventIndex_(orderedBoundaries.size()),
      operators_(nullptr) 
{
  if constexpr (!Master_matrix::Option_list::is_z2){
    operators_ = &(colSettings->operators);
  }

  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.reserve(orderedBoundaries.size());
  } else {
    pivotToColumnIndex_.resize(orderedBoundaries.size(), -1);
  }

  _initialize_U();
  _reduce();
}

template <class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(unsigned int numberOfColumns, 
                                           Column_settings* colSettings)
    : pair_opt(),
      swap_opt(),
      rep_opt(),
      reducedMatrixR_(numberOfColumns, colSettings),
      mirrorMatrixU_(numberOfColumns, colSettings),
      nextEventIndex_(0),
      operators_(nullptr) 
{
  if constexpr (!Master_matrix::Option_list::is_z2){
    operators_ = &(colSettings->operators);
  }

  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.reserve(numberOfColumns);
  } else {
    pivotToColumnIndex_.resize(numberOfColumns, -1);
  }
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _indexToBar().reserve(numberOfColumns);
  }
  if constexpr (Master_matrix::Option_list::has_vine_update) {
    swap_opt::_positionToRowIdx().reserve(numberOfColumns);
  }
}

template <class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(const RU_matrix& matrixToCopy, 
                                           Column_settings* colSettings)
    : pair_opt(static_cast<const pair_opt&>(matrixToCopy)),
      swap_opt(static_cast<const swap_opt&>(matrixToCopy)),
      rep_opt(static_cast<const rep_opt&>(matrixToCopy)),
      reducedMatrixR_(matrixToCopy.reducedMatrixR_, colSettings),
      mirrorMatrixU_(matrixToCopy.mirrorMatrixU_, colSettings),
      pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
      nextEventIndex_(matrixToCopy.nextEventIndex_),
      operators_(colSettings == nullptr ? matrixToCopy.operators_ : nullptr) 
{
  if constexpr (!Master_matrix::Option_list::is_z2){
    if (colSettings != nullptr) operators_ = &(colSettings->operators);
  }
}

template <class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(RU_matrix&& other) noexcept
    : pair_opt(std::move(static_cast<pair_opt&>(other))),
      swap_opt(std::move(static_cast<swap_opt&>(other))),
      rep_opt(std::move(static_cast<rep_opt&>(other))),
      reducedMatrixR_(std::move(other.reducedMatrixR_)),
      mirrorMatrixU_(std::move(other.mirrorMatrixU_)),
      pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
      nextEventIndex_(std::exchange(other.nextEventIndex_, 0)),
      operators_(std::exchange(other.operators_, nullptr)) 
{}

template <class Master_matrix>
template <class Boundary_type>
inline void RU_matrix<Master_matrix>::insert_boundary(const Boundary_type& boundary, dimension_type dim) 
{
  _insert_boundary(reducedMatrixR_.insert_boundary(boundary, dim));
}

template <class Master_matrix>
template <class Boundary_type>
inline void RU_matrix<Master_matrix>::insert_boundary(id_index faceIndex, 
                                                      const Boundary_type& boundary,
                                                      dimension_type dim) 
{
  //maps for possible shifting between column content and position indices used for birth events
  if constexpr (Master_matrix::Option_list::has_column_pairings && !Master_matrix::Option_list::has_vine_update){
    if (faceIndex != nextEventIndex_){
      pair_opt::idToPosition_.emplace(faceIndex, nextEventIndex_);
      if constexpr (Master_matrix::Option_list::has_removable_columns){
        pair_opt::RUM::map_.emplace(nextEventIndex_, faceIndex);
      }
    }
  }
  if constexpr (Master_matrix::Option_list::has_vine_update) {
    if (faceIndex != nextEventIndex_){
      swap_opt::_positionToRowIdx().emplace(nextEventIndex_, faceIndex);
      if (Master_matrix::Option_list::has_column_pairings){
        swap_opt::template RU_pairing<Master_matrix>::idToPosition_.emplace(faceIndex, nextEventIndex_);
      }
    }
  }
  _insert_boundary(reducedMatrixR_.insert_boundary(faceIndex, boundary, dim));
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Column_type& RU_matrix<Master_matrix>::get_column(index columnIndex,
                                                                                            bool inR) 
{
  if (inR) {
    return reducedMatrixR_.get_column(columnIndex);
  }
  return mirrorMatrixU_.get_column(columnIndex);
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::Row_type& RU_matrix<Master_matrix>::get_row(index rowIndex, bool inR) 
{
  static_assert(Master_matrix::Option_list::has_row_access, "'get_row' is not implemented for the chosen options.");

  if (inR) {
    return reducedMatrixR_.get_row(rowIndex);
  }
  return mirrorMatrixU_.get_row(rowIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::erase_empty_row(index rowIndex) 
{
  reducedMatrixR_.erase_empty_row(rowIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::remove_maximal_face(index columnIndex) 
{
  static_assert(Master_matrix::Option_list::has_removable_columns && Master_matrix::Option_list::has_vine_update,
                "'remove_maximal_face' is not implemented for the chosen options.");

  // TODO: is there an easy test to verify maximality even without row access?

  for (index curr = columnIndex; curr < nextEventIndex_ - 1; ++curr) {
    swap_opt::vine_swap(curr);
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
    id_index lastPivot = reducedMatrixR_.remove_last();
    if (lastPivot != static_cast<id_index>(-1)) pivotToColumnIndex_[lastPivot] = -1;
  }

  // if has_vine_update and has_column_pairings are both true,
  // then the element is already removed in _remove_last_in_barcode
  if constexpr (Master_matrix::Option_list::has_vine_update && !Master_matrix::Option_list::has_column_pairings) {
    swap_opt::_positionToRowIdx().erase(nextEventIndex_);
  }
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::dimension_type RU_matrix<Master_matrix>::get_max_dimension() const 
{
  return reducedMatrixR_.get_max_dimension();
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::index RU_matrix<Master_matrix>::get_number_of_columns() const 
{
  return reducedMatrixR_.get_number_of_columns();
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::dimension_type RU_matrix<Master_matrix>::get_column_dimension(
    index columnIndex) const 
{
  return reducedMatrixR_.get_column_dimension(columnIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex) 
{
  reducedMatrixR_.add_to(sourceColumnIndex, targetColumnIndex);
  //U transposed to avoid row operations
  if constexpr (Master_matrix::Option_list::has_vine_update)
    mirrorMatrixU_.add_to(targetColumnIndex, sourceColumnIndex);
  else
    mirrorMatrixU_.add_to(sourceColumnIndex, targetColumnIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::multiply_target_and_add_to(index sourceColumnIndex,
                                                                 const Field_element_type& coefficient,
                                                                 index targetColumnIndex) 
{
  reducedMatrixR_.multiply_target_and_add_to(sourceColumnIndex, coefficient, targetColumnIndex);
  mirrorMatrixU_.multiply_target_and_add_to(sourceColumnIndex, coefficient, targetColumnIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::multiply_source_and_add_to(const Field_element_type& coefficient,
                                                                 index sourceColumnIndex, 
                                                                 index targetColumnIndex) 
{
  reducedMatrixR_.multiply_source_and_add_to(coefficient, sourceColumnIndex, targetColumnIndex);
  mirrorMatrixU_.multiply_source_and_add_to(coefficient, sourceColumnIndex, targetColumnIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::zero_cell(index columnIndex, index rowIndex, bool inR) 
{
  if (inR) {
    return reducedMatrixR_.zero_cell(columnIndex, rowIndex);
  }
  return mirrorMatrixU_.zero_cell(columnIndex, rowIndex);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::zero_column(index columnIndex, bool inR) 
{
  if (inR) {
    return reducedMatrixR_.zero_column(columnIndex);
  }
  return mirrorMatrixU_.zero_column(columnIndex);
}

template <class Master_matrix>
inline bool RU_matrix<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex, bool inR) const 
{
  if (inR) {
    return reducedMatrixR_.is_zero_cell(columnIndex, rowIndex);
  }
  return mirrorMatrixU_.is_zero_cell(columnIndex, rowIndex);
}

template <class Master_matrix>
inline bool RU_matrix<Master_matrix>::is_zero_column(index columnIndex, bool inR) 
{
  if (inR) {
    return reducedMatrixR_.is_zero_column(columnIndex);
  }
  return mirrorMatrixU_.is_zero_column(columnIndex);
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::index RU_matrix<Master_matrix>::get_column_with_pivot(
    index faceIndex) const 
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return pivotToColumnIndex_.at(faceIndex);
  } else {
    return pivotToColumnIndex_[faceIndex];
  }
}

template <class Master_matrix>
inline typename RU_matrix<Master_matrix>::index RU_matrix<Master_matrix>::get_pivot(index columnIndex) 
{
  return reducedMatrixR_.get_column(columnIndex).get_pivot();
}

template <class Master_matrix>
inline RU_matrix<Master_matrix>& RU_matrix<Master_matrix>::operator=(const RU_matrix& other) 
{
  swap_opt::operator=(other);
  pair_opt::operator=(other);
  rep_opt::operator=(other);
  reducedMatrixR_ = other.reducedMatrixR_;
  mirrorMatrixU_ = other.mirrorMatrixU_;
  pivotToColumnIndex_ = other.pivotToColumnIndex_;
  nextEventIndex_ = other.nextEventIndex_;
  operators_ = other.operators_;
  return *this;
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
inline void RU_matrix<Master_matrix>::_insert_boundary(index currentIndex) 
{
  if constexpr (Master_matrix::Option_list::is_z2) {
    mirrorMatrixU_.insert_column({currentIndex});
  } else {
    mirrorMatrixU_.insert_column({{currentIndex, 1}});
  }

  if constexpr (!Master_matrix::Option_list::has_map_column_container) {
    id_index pivot = reducedMatrixR_.get_column(currentIndex).get_pivot();
    if (pivot != static_cast<id_index>(-1) && pivotToColumnIndex_.size() <= pivot)
      pivotToColumnIndex_.resize((pivot + 1) * 2, -1);
  }

  _reduce_last_column(currentIndex);
  ++nextEventIndex_;
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_initialize_U() 
{
  typename std::conditional<Master_matrix::Option_list::is_z2, index, std::pair<index, Field_element_type> >::type id;
  if constexpr (!Master_matrix::Option_list::is_z2) id.second = 1;

  for (id_index i = 0; i < reducedMatrixR_.get_number_of_columns(); i++) {
    if constexpr (Master_matrix::Option_list::is_z2)
      id = i;
    else
      id.first = i;
    mirrorMatrixU_.insert_column({id});
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce() 
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _indexToBar().reserve(reducedMatrixR_.get_number_of_columns());
  }

  for (index i = 0; i < reducedMatrixR_.get_number_of_columns(); i++) {
    if (!(reducedMatrixR_.is_zero_column(i))) {
      _reduce_column(i, i);
    } else {
      _add_bar(get_column_dimension(i), i);
    }
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce_last_column(index lastIndex) 
{
  if (reducedMatrixR_.get_column(lastIndex).is_empty()) {
    _add_bar(get_column_dimension(lastIndex), nextEventIndex_);
    return;
  }

  _reduce_column(lastIndex, nextEventIndex_);
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce_column(index target, index eventIndex)
{
  auto get_column_with_pivot_ = [&](id_index pivot) -> index {
    if (pivot == static_cast<id_index>(-1)) return -1;
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      auto it = pivotToColumnIndex_.find(pivot);
      if (it == pivotToColumnIndex_.end())
        return -1;
      else
        return it->second;
    } else {
      return pivotToColumnIndex_[pivot];
    }
  };

  Column_type& curr = reducedMatrixR_.get_column(target);
  id_index pivot = curr.get_pivot();
  index currIndex = get_column_with_pivot_(pivot);

  while (pivot != static_cast<id_index>(-1) && currIndex != static_cast<index>(-1)) {
    _reduce_column_by(target, currIndex);
    pivot = curr.get_pivot();
    currIndex = get_column_with_pivot_(pivot);
  }

  if (pivot != static_cast<id_index>(-1)) {
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
inline void RU_matrix<Master_matrix>::_reduce_column_by(index target, index source)
{
  Column_type& curr = reducedMatrixR_.get_column(target);
  if constexpr (Master_matrix::Option_list::is_z2) {
    curr += reducedMatrixR_.get_column(source);
    //to avoid having to do line operations during vineyards, U is transposed
    //TODO: explain this somewhere in the documentation...
    if constexpr (Master_matrix::Option_list::has_vine_update)
      mirrorMatrixU_.get_column(source) += mirrorMatrixU_.get_column(target);
    else
      mirrorMatrixU_.get_column(target) += mirrorMatrixU_.get_column(source);
  } else {
    Column_type& toadd = reducedMatrixR_.get_column(source);
    Field_element_type coef = toadd.get_pivot_value();
    coef = operators_->get_inverse(coef);
    operators_->multiply_inplace(coef, operators_->get_characteristic() - curr.get_pivot_value());

    curr.multiply_source_and_add(toadd, coef);
    mirrorMatrixU_.multiply_source_and_add_to(coef, source, target);
    // mirrorMatrixU_.get_column(target).multiply_source_and_add(mirrorMatrixU_.get_column(source), coef);
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_update_barcode(id_index birthPivot, pos_index death)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    if constexpr (Master_matrix::Option_list::has_vine_update)
      swap_opt::template RU_pairing<Master_matrix>::_update_barcode(birthPivot, death);
    else
      pair_opt::_update_barcode(birthPivot, death);
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_add_bar(dimension_type dim, pos_index birth) 
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    if constexpr (Master_matrix::Option_list::has_vine_update)
      swap_opt::template RU_pairing<Master_matrix>::_add_bar(dim, birth);
    else
      pair_opt::_add_bar(dim, birth);
  }
}

template <class Master_matrix>
inline void RU_matrix<Master_matrix>::_remove_last_in_barcode(pos_index eventIndex) 
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    if constexpr (Master_matrix::Option_list::has_vine_update)
      swap_opt::template RU_pairing<Master_matrix>::_remove_last(eventIndex);
    else
      pair_opt::_remove_last(eventIndex);
  }
}

template <class Master_matrix>
inline constexpr typename RU_matrix<Master_matrix>::bar_dictionary_type& RU_matrix<Master_matrix>::_indexToBar() 
{
  if constexpr (Master_matrix::Option_list::has_vine_update)
    return swap_opt::template RU_pairing<Master_matrix>::indexToBar_;
  else
    return pair_opt::indexToBar_;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_RU_MATRIX_H
