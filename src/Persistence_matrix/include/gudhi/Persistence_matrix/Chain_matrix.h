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
 * @file Chain_matrix.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Chain_matrix class.
 */

#ifndef PM_CHAIN_MATRIX_H
#define PM_CHAIN_MATRIX_H

#include <iostream>   //print() only
#include <set>
#include <map>
#include <stdexcept>
#include <vector>
#include <utility>    //std::swap, std::move & std::exchange
#include <algorithm>  //std::sort

#include <gudhi/Persistence_matrix/Id_to_index_overlay.h> //friend

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Chain_matrix Chain_matrix.h gudhi/Persistence_matrix/Chain_matrix.h
 * @ingroup persistence_matrix
 *
 * @brief %Matrix structure storing a compatible base of a filtered chain complex. See @cite zigzag.
 * The base is constructed from the boundaries of the cells in the complex. Allows the persistent homology to be
 * computed, as well as representative cycles. Supports vineyards (see @cite vineyards) and the removal 
 * of maximal cells while maintaining a valid barcode. Provides an access to its columns and rows.
 * 
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Chain_matrix : public Master_matrix::Matrix_dimension_option,
                     public Master_matrix::Chain_pairing_option,
                     public Master_matrix::Chain_vine_swap_option,
                     public Master_matrix::Chain_representative_cycles_option,
                     public Master_matrix::Matrix_row_access_option 
{
 public:
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element = typename Master_matrix::Element;                    /**< Type of an field element. */
  using Column = typename Master_matrix::Column;                            /**< Column type. */
  using Row = typename Master_matrix::Row;                                  /**< Row type, only necessary with row
                                                                                 access option. */
  using Entry = typename Master_matrix::Matrix_entry;                       /**< @ref Entry "Matrix entry" type. */
  using Entry_constructor = typename Master_matrix::Entry_constructor;      /**< Factory of @ref Entry classes. */
  using Column_settings = typename Master_matrix::Column_settings;          /**< Structure giving access to the columns
                                                                                 to necessary external classes. */
  using Boundary = typename Master_matrix::Boundary;                        /**< Type of an input column. */
  using Entry_representative = typename Master_matrix::Entry_representative;  /**< %Entry content representative. */
  using Index = typename Master_matrix::Index;                              /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;                        /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index;                      /**< @ref PosIdx index type. */
  using Dimension = typename Master_matrix::Dimension;                      /**< Dimension value type. */

  /**
   * @brief Constructs an empty matrix. Only available if @ref PersistenceMatrixOptions::has_column_pairings is
   * true or @ref PersistenceMatrixOptions::has_vine_update is false. Otherwise, birth and death comparators have
   * to be provided.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Chain_matrix(Column_settings* colSettings);
  /**
   * @brief Constructs a new matrix from the given ranges of @ref Matrix::Entry_representative. Each range corresponds
   * to a column  (the order of the ranges are preserved). The content of the ranges is assumed to be sorted by
   * increasing IDs. The IDs of the simplices are also assumed to be consecutive, ordered by filtration value, starting
   * with 0. Only available if @ref PersistenceMatrixOptions::has_column_pairings is true or
   * @ref PersistenceMatrixOptions::has_vine_update is false. Otherwise, birth and death
   * comparators have to be provided.
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
   * @ref Chain_matrix(unsigned int numberOfColumns, Column_settings* colSettings)
   * constructor preferably).
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  template <class Boundary_range = Boundary>
  Chain_matrix(const std::vector<Boundary_range>& orderedBoundaries, 
               Column_settings* colSettings);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns. Only available
   * if @ref PersistenceMatrixOptions::has_column_pairings is true or @ref PersistenceMatrixOptions::has_vine_update
   * is false. Otherwise, birth and death comparators have to be provided.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Chain_matrix(unsigned int numberOfColumns, Column_settings* colSettings);
  /**
   * @brief Constructs an empty matrix and stores the given comparators.
   *
   * @warning If @ref PersistenceMatrixOptions::has_vine_update is false, the comparators are not used.
   * And if @ref PersistenceMatrixOptions::has_vine_update is true, but
   * @ref PersistenceMatrixOptions::has_column_pairings is also true, the comparators are ignored and
   * the current barcode is used to compare birth and deaths. Therefore it is useless to provide them in those cases.
   * 
   * @tparam BirthComparatorFunction Type of the birth comparator: (@ref Pos_index, @ref Pos_index) -> bool
   * @tparam DeathComparatorFunction Type of the death comparator: (@ref Pos_index, @ref Pos_index) -> bool
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   * @param birthComparator Method taking two @ref PosIdx indices as input and returning true if and only if
   * the birth associated to the first position is strictly less than birth associated to
   * the second one with respect to some self defined order. It is used while swapping two unpaired or
   * two negative columns.
   * @param deathComparator Method taking two @ref PosIdx indices as input and returning true if and only if
   * the death associated to the first position is strictly less than death associated to
   * the second one with respect to some self defined order. It is used while swapping two positive but paired
   * columns.
   */
  template <typename BirthComparatorFunction, typename DeathComparatorFunction>
  Chain_matrix(Column_settings* colSettings, 
               const BirthComparatorFunction& birthComparator,
               const DeathComparatorFunction& deathComparator);
  /**
   * @brief Constructs a new matrix from the given ranges of @ref Matrix::Entry_representative. Each range corresponds
   * to a column (the order of the ranges are preserved). The content of the ranges is assumed to be sorted by
   * increasing IDs. The IDs of the simplices are also assumed to be consecutive, ordered by filtration value, starting
   * with 0.
   *
   * @warning If @ref PersistenceMatrixOptions::has_vine_update is false, the comparators are not used.
   * And if @ref PersistenceMatrixOptions::has_vine_update is true, but
   * @ref PersistenceMatrixOptions::has_column_pairings is also true, the comparators are ignored and
   * the current barcode is used to compare birth and deaths. Therefore it is useless to provide them in those cases.
   * 
   * @tparam BirthComparatorFunction Type of the birth comparator: (@ref Pos_index, @ref Pos_index) -> bool
   * @tparam DeathComparatorFunction Type of the death comparator: (@ref Pos_index, @ref Pos_index) -> bool
   * @tparam Boundary_range  Range type for @ref Matrix::Entry_representative ranges.
   * Assumed to have a begin(), end() and size() method.
   * @param orderedBoundaries Range of boundaries: @p orderedBoundaries is interpreted as a boundary matrix of a 
   * filtered **simplicial** complex, whose boundaries are ordered by filtration order. 
   * Therefore, `orderedBoundaries[i]` should store the boundary of the \f$ i^{th} \f$ simplex in the filtration,
   * as an ordered list of indices of its facets (again those indices correspond to their respective position
   * in the matrix). That is why the indices of the simplices are assumed to be consecutive and starting with 0 
   * (an empty boundary is interpreted as a vertex boundary and not as a non existing simplex). 
   * All dimensions up to the maximal dimension of interest have to be present. If only a higher dimension is of 
   * interest and not everything should be stored, then use the @ref insert_boundary method instead
   * (after creating the matrix with the @ref Chain_matrix(unsigned int, Column_settings*,
   * const BirthComparatorFunction&, const DeathComparatorFunction&) constructor preferably).
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   * @param birthComparator Method taking two @ref PosIdx indices as input and returning true if and only if
   * the birth associated to the first position is strictly less than birth associated to
   * the second one with respect to some self defined order. It is used while swapping two unpaired or
   * two negative columns.
   * @param deathComparator Method taking two @ref PosIdx indices as input and returning true if and only if
   * the death associated to the first position is strictly less than death associated to
   * the second one with respect to some self defined order. It is used while swapping two positive but paired
   * columns.
   */
  template <typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_range = Boundary>
  Chain_matrix(const std::vector<Boundary_range>& orderedBoundaries, 
               Column_settings* colSettings, 
               const BirthComparatorFunction& birthComparator,
               const DeathComparatorFunction& deathComparator);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   *
   * @warning If @ref PersistenceMatrixOptions::has_vine_update is false, the comparators are not used.
   * And if @ref PersistenceMatrixOptions::has_vine_update is true, but
   * @ref PersistenceMatrixOptions::has_column_pairings is also true, the comparators are ignored and
   * the current barcode is used to compare birth and deaths. Therefore it is useless to provide them in those cases.
   * 
   * @tparam BirthComparatorFunction Type of the birth comparator: (@ref Pos_index, @ref Pos_index) -> bool
   * @tparam DeathComparatorFunction Type of the death comparator: (@ref Pos_index, @ref Pos_index) -> bool
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   * @param birthComparator Method taking two @ref PosIdx indices as input and returning true if and only if
   * the birth associated to the first position is strictly less than birth associated to
   * the second one with respect to some self defined order. It is used while swapping two unpaired or
   * two negative columns.
   * @param deathComparator Method taking two @ref PosIdx indices as input and returning true if and only if
   * the death associated to the first position is strictly less than death associated to
   * the second one with respect to some self defined order. It is used while swapping two positive but paired
   * columns.
   */
  template <typename BirthComparatorFunction, typename DeathComparatorFunction>
  Chain_matrix(unsigned int numberOfColumns, 
               Column_settings* colSettings,
               const BirthComparatorFunction& birthComparator, 
               const DeathComparatorFunction& deathComparator);
  /**
   * @brief Copy constructor. If @p colSettings is not a null pointer, its value is kept
   * instead of the one in the copied matrix.
   * 
   * @param matrixToCopy Matrix to copy.
   * @param colSettings Either a pointer to an existing setting structure for the columns or a null pointer.
   * The structure should contain all the necessary external classes specifically necessary for the choosen column type,
   * such as custom allocators. If null pointer, the pointer stored in @p matrixToCopy is used instead.
   */
  Chain_matrix(const Chain_matrix& matrixToCopy, 
               Column_settings* colSettings = nullptr);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Chain_matrix(Chain_matrix&& other) noexcept;

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
   * When inserted, the given boundary is reduced and from the reduction process, the column is deduced in the form of:
   * `IDIdx + linear combination of older column IDIdxs`. If the barcode is stored, it will be updated.
   * 
   * @tparam Boundary_range Range of @ref Matrix::Entry_representative. Assumed to have a begin(), end() and size()
   * method.
   * @param boundary Boundary generating the new column. The content should be ordered by ID.
   * @param dim Dimension of the cell whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return The @ref MatIdx indices of the unpaired chains used to reduce the boundary.
   */
  template <class Boundary_range = Boundary>
  std::vector<Entry_representative> insert_boundary(const Boundary_range& boundary, Dimension dim = -1);
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
   * @param cellID @ref IDIdx index to use to identify the new cell.
   * @param boundary Boundary generating the new column. The indices of the boundary have to correspond to the 
   * @p cellID values of precedent calls of the method for the corresponding cells and should be ordered in 
   * increasing order.
   * @param dim Dimension of the cell whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return The @ref MatIdx index of the inserted boundary.
   */
  template <class Boundary_range = Boundary>
  std::vector<Entry_representative> insert_boundary(ID_index cellID,
                                                    const Boundary_range& boundary,
                                                    Dimension dim = -1);
  /**
   * @brief Returns the column at the given @ref MatIdx index.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param columnIndex @ref MatIdx index of the column to return.
   * @return Reference to the column.
   */
  Column& get_column(Index columnIndex);
  /**
   * @brief Returns the column at the given @ref MatIdx index.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param columnIndex @ref MatIdx index of the column to return.
   * @return Const reference to the column.
   */
  const Column& get_column(Index columnIndex) const;
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns and
   * @ref PersistenceMatrixOptions::has_vine_update are true, as well as,
   * @ref PersistenceMatrixOptions::has_map_column_container and @ref PersistenceMatrixOptions::has_column_pairings.
   * Assumes that the cell is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., the matrix is still a compatible bases of the chain complex in the sense of @cite zigzag).
   * The maximality of the cell is not verified.
   * Also updates the barcode if it is stored.
   *
   * Note that using the other version of the method could perform better depending on how the data is 
   * maintained on the side of the user, that is, if providing the second parameter is easy.
   *
   * See also @ref remove_last.
   * 
   * @param cellID @ref IDIdx index of the cell to remove
   */
  void remove_maximal_cell(ID_index cellID);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns,
   * @ref PersistenceMatrixOptions::has_vine_update and @ref PersistenceMatrixOptions::has_map_column_container
   * are true.
   * Assumes that the cell is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., it is still a compatible bases of the chain complex in the sense of @cite zigzag).
   * The maximality of the cell is not verified.
   * Also updates the barcode if it is stored.
   *
   * To maintain the compatibility, vine swaps are done to move the cell up to the end of the filtration. Once at 
   * the end, the removal is trivial. But for @ref chainmatrix "chain matrices", swaps do not actually swap the position
   * of the column every time, so the cells appearing after @p cellID in the filtration have to be searched first within
   * the matrix. If the user has an easy access to the @ref IDIdx of the cells in the order of filtration, passing them
   * by argument with @p columnsToSwap allows to skip a linear search process. Typically, if the user knows that the
   * cell he wants to remove is already the last cell of the filtration, calling
   * @ref remove_maximal_cell(ID_index cellIndex, const std::vector<ID_index>& columnsToSwap)
   * "remove_maximal_cell(cellID, {})" will be faster than @ref remove_last().
   *
   * See also @ref remove_last.
   * 
   * @param cellID @ref IDIdx index of the cell to remove
   * @param columnsToSwap Vector of @ref IDIdx indices of the cells coming after @p cellID in the filtration.
   */
  void remove_maximal_cell(ID_index cellID, const std::vector<ID_index>& columnsToSwap);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns is true and,
   * if @ref PersistenceMatrixOptions::has_map_column_container is true or
   * @ref PersistenceMatrixOptions::has_vine_update is false.
   * Removes the last cell in the filtration from the matrix and updates the barcode if it is stored.
   *
   * See also @ref remove_maximal_cell.
   *
   * @warning If @ref PersistenceMatrixOptions::has_vine_update is true, the last cell does not have to
   * be at the end of the matrix container and therefore has to be searched first. In this case, if the user
   * already knows the @ref IDIdx of the last cell, calling
   * @ref remove_maximal_cell(ID_index cellIndex, const std::vector<ID_index>& columnsToSwap)
   * "remove_maximal_cell(cellID, {})" instead allows to skip the search.
   */
  void remove_last();

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
   * @ref chainmatrix "chain matrix". For example, a right-to-left addition could corrupt the computation
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
   * @ref chainmatrix "chain matrix". For example, a right-to-left addition could corrupt the computation
   * of the barcode if done blindly. So should be used with care.
   * 
   * @param sourceColumnIndex @ref MatIdx index of the source column.
   * @param coefficient Value to multiply.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  void multiply_target_and_add_to(Index sourceColumnIndex, 
                                  const Field_element& coefficient,
                                  Index targetColumnIndex);
  /**
   * @brief Multiplies the source column with the coefficient before adding it to the target column.
   * That is: `targetColumn += (coefficient * sourceColumn)`. The source column will **not** be modified.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of a
   * @ref chainmatrix "chain matrix". For example, a right-to-left addition could corrupt the computation
   * of the barcode if done blindly. So should be used with care.
   * 
   * @param coefficient Value to multiply.
   * @param sourceColumnIndex @ref MatIdx index of the source column.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  void multiply_source_and_add_to(const Field_element& coefficient, 
                                  Index sourceColumnIndex,
                                  Index targetColumnIndex);

  /**
   * @brief Indicates if the entry at given coordinates has value zero.
   * 
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @return true If the entry has value zero.
   * @return false Otherwise.
   */
  bool is_zero_entry(Index columnIndex, ID_index rowIndex) const;
  /**
   * @brief Indicates if the column at given index has value zero. Note that if the matrix is valid, this method
   * should always return false.
   * 
   * @param columnIndex @ref MatIdx index of the column.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(Index columnIndex);

  /**
   * @brief Returns the column with given @ref rowindex "row index" as pivot. Assumes that the pivot exists.
   * 
   * @param cellID @ref rowindex "Row index" of the pivot.
   * @return @ref MatIdx index of the column with the given pivot.
   */
  Index get_column_with_pivot(ID_index cellID) const;
  /**
   * @brief Returns the @ref rowindex "row index" of the pivot of the given column.
   * 
   * @param columnIndex @ref MatIdx index of the column
   * @return The @ref rowindex "row index" of the pivot.
   */
  ID_index get_pivot(Index columnIndex);

  /**
   * @brief Resets the matrix to an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  void reset(Column_settings* colSettings) {
    matrix_.clear();
    pivotToColumnIndex_.clear();
    nextIndex_ = 0;
    colSettings_ = colSettings;
  }

  /**
   * @brief Assign operator.
   */
  Chain_matrix& operator=(const Chain_matrix& other);
  /**
   * @brief Swap operator.
   */
  friend void swap(Chain_matrix& matrix1, Chain_matrix& matrix2) {
    swap(static_cast<typename Master_matrix::Matrix_dimension_option&>(matrix1),
         static_cast<typename Master_matrix::Matrix_dimension_option&>(matrix2));
    swap(static_cast<typename Master_matrix::Chain_pairing_option&>(matrix1),
         static_cast<typename Master_matrix::Chain_pairing_option&>(matrix2));
    swap(static_cast<typename Master_matrix::Chain_vine_swap_option&>(matrix1),
         static_cast<typename Master_matrix::Chain_vine_swap_option&>(matrix2));
    swap(static_cast<typename Master_matrix::Chain_representative_cycles_option&>(matrix1),
         static_cast<typename Master_matrix::Chain_representative_cycles_option&>(matrix2));
    matrix1.matrix_.swap(matrix2.matrix_);
    matrix1.pivotToColumnIndex_.swap(matrix2.pivotToColumnIndex_);
    std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
    std::swap(matrix1.colSettings_, matrix2.colSettings_);

    if constexpr (Master_matrix::Option_list::has_row_access) {
      swap(static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix1),
           static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix2));
    }
  }

  void print() const;  // for debug

  friend class Id_to_index_overlay<Chain_matrix<Master_matrix>, Master_matrix>;

 private:
  using Dim_opt = typename Master_matrix::Matrix_dimension_option;
  using Swap_opt = typename Master_matrix::Chain_vine_swap_option;
  using Pair_opt = typename Master_matrix::Chain_pairing_option;
  using Rep_opt = typename Master_matrix::Chain_representative_cycles_option;
  using RA_opt = typename Master_matrix::Matrix_row_access_option;
  using Column_container = typename Master_matrix::Column_container;
  using Dictionary = typename Master_matrix::template Dictionary<Index>;
  using Barcode = typename Master_matrix::Barcode;
  using Bar_dictionary = typename Master_matrix::Bar_dictionary;
  using Tmp_column = typename std::conditional<Master_matrix::Option_list::is_z2,
                                               std::set<ID_index>,
                                               std::map<ID_index, Field_element>
                                              >::type;

  Column_container matrix_;       /**< Column container. */
  Dictionary pivotToColumnIndex_; /**< Map from @ref IDIdx to @ref MatIdx index. */
  Index nextIndex_;               /**< Next unused column index. */
  Column_settings* colSettings_;  /**< Entry factory. */

  template <class Boundary_range>
  std::vector<Entry_representative> _reduce_boundary(ID_index cellID, const Boundary_range& boundary, Dimension dim);
  void _reduce_by_G(Tmp_column& column, std::vector<Entry_representative>& chainsInH, Index currentPivot);
  void _reduce_by_F(Tmp_column& column, std::vector<Entry_representative>& chainsInF, Index currentPivot);
  void _build_from_H(ID_index cellID, Tmp_column& column, std::vector<Entry_representative>& chainsInH);
  void _update_largest_death_in_F(const std::vector<Entry_representative>& chainsInF);
  void _insert_chain(const Tmp_column& column, Dimension dimension);
  void _insert_chain(const Tmp_column& column, Dimension dimension, Index pair);
  void _add_to(const Column& column, Tmp_column& set, unsigned int coef);
  template <typename F>
  void _add_to(Column& target, F&& addition);
  void _remove_last(Index lastIndex);
  void _update_barcode(Pos_index birth);
  void _add_bar(Dimension dim);
  template <class Container>
  void _container_insert(const Container& column, Index pos, Dimension dim);
  void _container_insert(const Column& column, [[maybe_unused]] Index pos = 0);

  constexpr Barcode& _barcode();
  constexpr Bar_dictionary& _indexToBar();
  constexpr Pos_index& _nextPosition();
};

template <class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(Column_settings* colSettings)
    : Dim_opt(-1),
      Pair_opt(),
      Swap_opt(),
      Rep_opt(),
      RA_opt(),
      nextIndex_(0),
      colSettings_(colSettings)
{}

template <class Master_matrix>
template <class Boundary_range>
inline Chain_matrix<Master_matrix>::Chain_matrix(const std::vector<Boundary_range>& orderedBoundaries,
                                                 Column_settings* colSettings)
    : Dim_opt(-1),
      Pair_opt(),
      Swap_opt(),
      Rep_opt(),
      RA_opt(orderedBoundaries.size()),
      nextIndex_(0),
      colSettings_(colSettings)
{
  matrix_.reserve(orderedBoundaries.size());
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.reserve(orderedBoundaries.size());
  } else {
    pivotToColumnIndex_.resize(orderedBoundaries.size(), -1);
  }

  for (const Boundary_range& b : orderedBoundaries) {
    insert_boundary(b);
  }
}

template <class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(unsigned int numberOfColumns, 
                                                 Column_settings* colSettings)
    : Dim_opt(-1),
      Pair_opt(),
      Swap_opt(),
      Rep_opt(),
      RA_opt(numberOfColumns),
      nextIndex_(0),
      colSettings_(colSettings)
{
  matrix_.reserve(numberOfColumns);
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.reserve(numberOfColumns);
  } else {
    pivotToColumnIndex_.resize(numberOfColumns, -1);
  }
}

template <class Master_matrix>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Chain_matrix<Master_matrix>::Chain_matrix(Column_settings* colSettings,
                                                 const BirthComparatorFunction& birthComparator,
                                                 const DeathComparatorFunction& deathComparator)
    : Dim_opt(-1),
      Pair_opt(),
      Swap_opt(birthComparator, deathComparator),
      Rep_opt(),
      RA_opt(),
      nextIndex_(0),
      colSettings_(colSettings)
{}

template <class Master_matrix>
template <typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_range>
inline Chain_matrix<Master_matrix>::Chain_matrix(const std::vector<Boundary_range>& orderedBoundaries,
                                                 Column_settings* colSettings,
                                                 const BirthComparatorFunction& birthComparator,
                                                 const DeathComparatorFunction& deathComparator)
    : Dim_opt(-1),
      Pair_opt(),
      Swap_opt(birthComparator, deathComparator),
      Rep_opt(),
      RA_opt(orderedBoundaries.size()),
      nextIndex_(0),
      colSettings_(colSettings)
{
  matrix_.reserve(orderedBoundaries.size());
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.reserve(orderedBoundaries.size());
  } else {
    pivotToColumnIndex_.resize(orderedBoundaries.size(), -1);
  }
  for (const Boundary_range& b : orderedBoundaries) {
    insert_boundary(b);
  }
}

template <class Master_matrix>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Chain_matrix<Master_matrix>::Chain_matrix(unsigned int numberOfColumns, 
                                                 Column_settings* colSettings,
                                                 const BirthComparatorFunction& birthComparator,
                                                 const DeathComparatorFunction& deathComparator)
    : Dim_opt(-1),
      Pair_opt(),
      Swap_opt(birthComparator, deathComparator),
      Rep_opt(),
      RA_opt(numberOfColumns),
      nextIndex_(0),
      colSettings_(colSettings)
{
  matrix_.reserve(numberOfColumns);
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.reserve(numberOfColumns);
  } else {
    pivotToColumnIndex_.resize(numberOfColumns, -1);
  }
}

template <class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(const Chain_matrix& matrixToCopy, Column_settings* colSettings)
    : Dim_opt(static_cast<const Dim_opt&>(matrixToCopy)),
      Pair_opt(static_cast<const Pair_opt&>(matrixToCopy)),
      Swap_opt(static_cast<const Swap_opt&>(matrixToCopy)),
      Rep_opt(static_cast<const Rep_opt&>(matrixToCopy)),
      RA_opt(static_cast<const RA_opt&>(matrixToCopy)),
      pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
      nextIndex_(matrixToCopy.nextIndex_),
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
inline Chain_matrix<Master_matrix>::Chain_matrix(Chain_matrix&& other) noexcept
    : Dim_opt(std::move(static_cast<Dim_opt&>(other))),
      Pair_opt(std::move(static_cast<Pair_opt&>(other))),
      Swap_opt(std::move(static_cast<Swap_opt&>(other))),
      Rep_opt(std::move(static_cast<Rep_opt&>(other))),
      RA_opt(std::move(static_cast<RA_opt&>(other))),
      matrix_(std::move(other.matrix_)),
      pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
      nextIndex_(std::exchange(other.nextIndex_, 0)),
      colSettings_(std::exchange(other.colSettings_, nullptr)) 
{}

template <class Master_matrix>
template <class Boundary_range>
inline std::vector<typename Master_matrix::Entry_representative> Chain_matrix<Master_matrix>::insert_boundary(
    const Boundary_range& boundary, Dimension dim) 
{
  return insert_boundary(nextIndex_, boundary, dim);
}

template <class Master_matrix>
template <class Boundary_range>
inline std::vector<typename Master_matrix::Entry_representative> Chain_matrix<Master_matrix>::insert_boundary(
    ID_index cellID, const Boundary_range& boundary, Dimension dim) 
{
  if constexpr (!Master_matrix::Option_list::has_map_column_container) {
    if (pivotToColumnIndex_.size() <= cellID) {
      pivotToColumnIndex_.resize(cellID * 2 + 1, -1);
    }
  }

  if constexpr (Master_matrix::Option_list::has_vine_update && Master_matrix::Option_list::has_column_pairings) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      Swap_opt::CP::pivotToPosition_.try_emplace(cellID, _nextPosition());
    } else {
      if (Swap_opt::CP::pivotToPosition_.size() <= cellID)
        Swap_opt::CP::pivotToPosition_.resize(pivotToColumnIndex_.size(), -1);
      Swap_opt::CP::pivotToPosition_[cellID] = _nextPosition();
    }
  }

  if constexpr (Master_matrix::Option_list::has_matrix_maximal_dimension_access) {
    Dim_opt::update_up(dim == static_cast<Dimension>(-1) ? (boundary.size() == 0 ? 0 : boundary.size() - 1) : dim);
  }

  return _reduce_boundary(cellID, boundary, dim);
}

template <class Master_matrix>
inline typename Chain_matrix<Master_matrix>::Column& Chain_matrix<Master_matrix>::get_column(Index columnIndex) 
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.at(columnIndex);
  } else {
    return matrix_[columnIndex];
  }
}

template <class Master_matrix>
inline const typename Chain_matrix<Master_matrix>::Column& Chain_matrix<Master_matrix>::get_column(
    Index columnIndex) const
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.at(columnIndex);
  } else {
    return matrix_[columnIndex];
  }
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::remove_maximal_cell(ID_index cellID)
{
  static_assert(Master_matrix::Option_list::has_removable_columns,
                "'remove_maximal_cell' is not implemented for the chosen options.");
  static_assert(Master_matrix::Option_list::has_map_column_container && 
                    Master_matrix::Option_list::has_vine_update &&
                    Master_matrix::Option_list::has_column_pairings,
                "'remove_maximal_cell' is not implemented for the chosen options.");

  // TODO: find simple test to verify that col at columnIndex is maximal even without row access.

  const auto& pivotToPosition = Swap_opt::CP::pivotToPosition_;
  auto it = pivotToPosition.find(cellID);
  if (it == pivotToPosition.end()) return;  // cell does not exists. TODO: put an assert instead?
  Pos_index startPos = it->second;
  Index startIndex = pivotToColumnIndex_.at(cellID);

  if (startPos != _nextPosition() - 1) {
    std::vector<Index> colToSwap;
    colToSwap.reserve(matrix_.size());

    for (auto& p : pivotToPosition) {
      if (p.second > startPos) colToSwap.push_back(pivotToColumnIndex_.at(p.first));
    }
    std::sort(colToSwap.begin(), colToSwap.end(), [&](Index c1, Index c2) {
      return pivotToPosition.at(get_pivot(c1)) < pivotToPosition.at(get_pivot(c2));
    });

    for (Index i : colToSwap) {
      startIndex = Swap_opt::vine_swap(startIndex, i);
    }
  }

  _remove_last(startIndex);
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::remove_maximal_cell(ID_index cellID,
                                                             const std::vector<ID_index>& columnsToSwap)
{
  static_assert(Master_matrix::Option_list::has_removable_columns,
                "'remove_maximal_cell' is not implemented for the chosen options.");
  static_assert(Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_vine_update,
                "'remove_maximal_cell' is not implemented for the chosen options.");

  // TODO: find simple test to verify that col at columnIndex is maximal even without row access.

  Index startIndex = pivotToColumnIndex_.at(cellID);

  for (ID_index i : columnsToSwap) {
    startIndex = Swap_opt::vine_swap(startIndex, pivotToColumnIndex_.at(i));
  }

  _remove_last(startIndex);
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::remove_last()
{
  static_assert(Master_matrix::Option_list::has_removable_columns,
                "'remove_last' is not implemented for the chosen options.");
  static_assert(Master_matrix::Option_list::has_map_column_container || !Master_matrix::Option_list::has_vine_update,
                "'remove_last' is not implemented for the chosen options.");

  if (nextIndex_ == 0 || matrix_.empty()) return;  // empty matrix

  if constexpr (Master_matrix::Option_list::has_vine_update) {
    // careful: linear because of the search of the last index. It is better to keep track of the @ref IDIdx index
    // of the last column while performing swaps (or the @ref MatIdx with the return values of `vine_swap` + get_pivot)
    // and then call `remove_maximal_cell` with it and an empty `columnsToSwap`.

    ID_index pivot = 0;
    Index colIndex = 0;
    for (auto& p : pivotToColumnIndex_) {
      if (p.first > pivot) {  // pivots have to be strictly increasing in order of filtration
        pivot = p.first;
        colIndex = p.second;
      }
    }
    _remove_last(colIndex);
  } else {
    _remove_last(nextIndex_ - 1);
  }
}

template <class Master_matrix>
inline typename Chain_matrix<Master_matrix>::Index Chain_matrix<Master_matrix>::get_number_of_columns() const
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return matrix_.size();
  } else {
    return nextIndex_;  // matrix could have been resized much bigger while insert
  }
}

template <class Master_matrix>
inline typename Chain_matrix<Master_matrix>::Dimension Chain_matrix<Master_matrix>::get_column_dimension(
    Index columnIndex) const 
{
  return get_column(columnIndex).get_dimension();
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::add_to(Index sourceColumnIndex, Index targetColumnIndex)
{
  auto& col = get_column(targetColumnIndex);
  _add_to(col, [&]() { col += get_column(sourceColumnIndex); });
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::multiply_target_and_add_to(Index sourceColumnIndex,
                                                                    const Field_element& coefficient,
                                                                    Index targetColumnIndex)
{
  auto& col = get_column(targetColumnIndex);
  _add_to(col, [&]() { col.multiply_target_and_add(coefficient, get_column(sourceColumnIndex)); });
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::multiply_source_and_add_to(const Field_element& coefficient,
                                                                    Index sourceColumnIndex,
                                                                    Index targetColumnIndex)
{
  auto& col = get_column(targetColumnIndex);
  _add_to(col, [&]() { col.multiply_source_and_add(get_column(sourceColumnIndex), coefficient); });
}

template <class Master_matrix>
inline bool Chain_matrix<Master_matrix>::is_zero_entry(Index columnIndex, ID_index rowIndex) const 
{
  return !get_column(columnIndex).is_non_zero(rowIndex);
}

template <class Master_matrix>
inline bool Chain_matrix<Master_matrix>::is_zero_column(Index columnIndex) 
{
  return get_column(columnIndex).is_empty();
}

template <class Master_matrix>
inline typename Chain_matrix<Master_matrix>::Index Chain_matrix<Master_matrix>::get_column_with_pivot(
    ID_index cellID) const
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return pivotToColumnIndex_.at(cellID);
  } else {
    return pivotToColumnIndex_[cellID];
  }
}

template <class Master_matrix>
inline typename Chain_matrix<Master_matrix>::ID_index Chain_matrix<Master_matrix>::get_pivot(Index columnIndex)
{
  return get_column(columnIndex).get_pivot();
}

template <class Master_matrix>
inline Chain_matrix<Master_matrix>& Chain_matrix<Master_matrix>::operator=(const Chain_matrix& other)
{
  Dim_opt::operator=(other);
  Swap_opt::operator=(other);
  Pair_opt::operator=(other);
  Rep_opt::operator=(other);
  matrix_.clear();
  pivotToColumnIndex_ = other.pivotToColumnIndex_;
  nextIndex_ = other.nextIndex_;
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
inline void Chain_matrix<Master_matrix>::print() const
{
  std::cout << "Column Matrix:\n";
  if constexpr (!Master_matrix::Option_list::has_map_column_container) {
    for (ID_index i = 0; i < pivotToColumnIndex_.size() && pivotToColumnIndex_[i] != static_cast<Index>(-1); ++i) {
      Index pos = pivotToColumnIndex_[i];
      const Column& col = matrix_[pos];
      for (const auto& entry : col) {
        std::cout << entry.get_row_index() << " ";
      }
      std::cout << "(" << i << ", " << pos << ")\n";
    }
    if constexpr (Master_matrix::Option_list::has_row_access) {
      std::cout << "\n";
      std::cout << "Row Matrix:\n";
      for (ID_index i = 0; i < pivotToColumnIndex_.size() && pivotToColumnIndex_[i] != static_cast<Index>(-1); ++i) {
        Index pos = pivotToColumnIndex_[i];
        const Row& row = RA_opt::get_row(pos);
        for (const auto& entry : row) {
          std::cout << entry.get_column_index() << " ";
        }
        std::cout << "(" << i << ", " << pos << ")\n";
      }
    }
  } else {
    for (const auto& p : pivotToColumnIndex_) {
      const Column& col = matrix_.at(p.second);
      for (const auto& entry : col) {
        std::cout << entry.get_row_index() << " ";
      }
      std::cout << "(" << p.first << ", " << p.second << ")\n";
    }
    if constexpr (Master_matrix::Option_list::has_row_access) {
      std::cout << "\n";
      std::cout << "Row Matrix:\n";
      for (const auto& p : pivotToColumnIndex_) {
        const Row& row = RA_opt::get_row(p.first);
        for (const auto& entry : row) {
          std::cout << entry.get_column_index() << " ";
        }
        std::cout << "(" << p.first << ", " << p.second << ")\n";
      }
    }
  }
  std::cout << "\n";
}

template <class Master_matrix>
template <class Boundary_range>
inline std::vector<typename Master_matrix::Entry_representative> Chain_matrix<Master_matrix>::_reduce_boundary(
    ID_index cellID, const Boundary_range& boundary, Dimension dim)
{
  Tmp_column column(boundary.begin(), boundary.end());
  if (dim == static_cast<Dimension>(-1)) dim = boundary.begin() == boundary.end() ? 0 : boundary.size() - 1;
  std::vector<Entry_representative> chainsInH;  // for corresponding indices in H (paired columns)
  std::vector<Entry_representative> chainsInF;  // for corresponding indices in F (unpaired, essential columns)

  auto get_last = [&column]() {
    if constexpr (Master_matrix::Option_list::is_z2)
      return *(column.rbegin());
    else
      return column.rbegin()->first;
  };

  if (boundary.begin() == boundary.end()) {
    if constexpr (Master_matrix::Option_list::is_z2)
      column.insert(cellID);
    else
      column.emplace(cellID, 1);
    _insert_chain(column, dim);
    return chainsInF;
  }

  Index currentIndex = get_column_with_pivot(get_last());

  while (get_column(currentIndex).is_paired()) {
    _reduce_by_G(column, chainsInH, currentIndex);

    if (column.empty()) {
      // produce the sum of all col_h in chains_in_H
      _build_from_H(cellID, column, chainsInH);
      // create a new cycle (in F) sigma - \sum col_h
      _insert_chain(column, dim);
      return chainsInF;
    }

    currentIndex = get_column_with_pivot(get_last());
  }

  while (!column.empty()) {
    currentIndex = get_column_with_pivot(get_last());

    if (!get_column(currentIndex).is_paired()) {
      // only fills currentEssentialCycleIndices if Z2 coefficients, so chainsInF remains empty
      _reduce_by_F(column, chainsInF, currentIndex);
    } else {
      _reduce_by_G(column, chainsInH, currentIndex);
    }
  }

  _update_largest_death_in_F(chainsInF);

  // Compute the new column cellID + \sum col_h, for col_h in chains_in_H
  _build_from_H(cellID, column, chainsInH);

  // Create and insert (\sum col_h) + sigma (in H, paired with chain_fp) in matrix_
  if constexpr (Master_matrix::Option_list::is_z2)
    _insert_chain(column, dim, chainsInF[0]);
  else
    _insert_chain(column, dim, chainsInF[0].first);

  return chainsInF;
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_reduce_by_G(Tmp_column& column,
                                                      std::vector<Entry_representative>& chainsInH,
                                                      Index currentIndex)
{
  Column& col = get_column(currentIndex);
  if constexpr (Master_matrix::Option_list::is_z2) {
    _add_to(col, column, 1u);                           // Reduce with the column col_g
    chainsInH.push_back(col.get_paired_chain_index());  // keep the col_h with which col_g is paired
  } else {
    Field_element coef = col.get_pivot_value();
    auto& operators = colSettings_->operators;
    coef = operators.get_inverse(coef);
    operators.multiply_inplace(coef, operators.get_characteristic() - column.rbegin()->second);

    _add_to(col, column, coef);                                       // Reduce with the column col_g
    chainsInH.emplace_back(col.get_paired_chain_index(), coef); // keep the col_h with which col_g is paired
  }
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_reduce_by_F(Tmp_column& column, 
                                                      std::vector<Entry_representative>& chainsInF,
                                                      Index currentIndex) 
{
  Column& col = get_column(currentIndex);
  if constexpr (Master_matrix::Option_list::is_z2) {
    _add_to(col, column, 1u);  // Reduce with the column col_g
    chainsInF.push_back(currentIndex);
  } else {
    Field_element coef = col.get_pivot_value();
    auto& operators = colSettings_->operators;
    coef = operators.get_inverse(coef);
    operators.multiply_inplace(coef, operators.get_characteristic() - column.rbegin()->second);

    _add_to(col, column, coef);  // Reduce with the column col_g
    chainsInF.emplace_back(currentIndex, operators.get_characteristic() - coef);
  }
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_build_from_H(ID_index cellID,
                                                       Tmp_column& column,
                                                       std::vector<Entry_representative>& chainsInH)
{
  if constexpr (Master_matrix::Option_list::is_z2) {
    column.insert(cellID);
    for (Index idx_h : chainsInH) {
      _add_to(get_column(idx_h), column, 1u);
    }
  } else {
    column.emplace(cellID, 1);
    for (std::pair<Index, Field_element>& idx_h : chainsInH) {
      _add_to(get_column(idx_h.first), column, idx_h.second);
    }
  }
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_update_largest_death_in_F(const std::vector<Entry_representative>& chainsInF) 
{
  if constexpr (Master_matrix::Option_list::is_z2) {
    Index toUpdate = chainsInF[0];
    for (auto other_col_it = chainsInF.begin() + 1; other_col_it != chainsInF.end(); ++other_col_it) {
      add_to(*other_col_it, toUpdate);
    }
  } else {
    Index toUpdate = chainsInF[0].first;
    get_column(toUpdate) *= chainsInF[0].second;
    for (auto other_col_it = chainsInF.begin() + 1; other_col_it != chainsInF.end(); ++other_col_it) {
      multiply_source_and_add_to(other_col_it->second, other_col_it->first, toUpdate);
    }
  }
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_insert_chain(const Tmp_column& column, Dimension dimension) 
{
  _container_insert(column, nextIndex_, dimension);
  _add_bar(dimension);

  ++nextIndex_;
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_insert_chain(const Tmp_column& column, Dimension dimension, Index pair)
{
  // true when no vine updates and if nextIndex_ is updated in remove_last for special case of no vines
  // because then @ref PosIdx == @ref MatIdx
  Pos_index pairPos = pair;

  _container_insert(column, nextIndex_, dimension);

  get_column(nextIndex_).assign_paired_chain(pair);
  auto& pairCol = get_column(pair);
  pairCol.assign_paired_chain(nextIndex_);

  if constexpr (Master_matrix::Option_list::has_column_pairings && Master_matrix::Option_list::has_vine_update) {
    pairPos = Swap_opt::CP::pivotToPosition_[pairCol.get_pivot()];
  }

  _update_barcode(pairPos);

  ++nextIndex_;
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_add_to(const Column& column,
                                                 Tmp_column& set,
                                                 [[maybe_unused]] unsigned int coef) 
{
  if constexpr (Master_matrix::Option_list::is_z2) {
    std::pair<typename std::set<Index>::iterator, bool> res_insert;
    for (const Entry& entry : column) {
      res_insert = set.insert(entry.get_row_index());
      if (!res_insert.second) {
        set.erase(res_insert.first);
      }
    }
  } else {
    auto& operators = colSettings_->operators;
    for (const Entry& entry : column) {
      auto res = set.emplace(entry.get_row_index(), entry.get_element());
      if (res.second){
        operators.multiply_inplace(res.first->second, coef);
      } else {
        operators.multiply_and_add_inplace_back(entry.get_element(), coef, res.first->second);
        if (res.first->second == Field_operators::get_additive_identity()) {
          set.erase(res.first);
        }
      }
    }
  }
}

template <class Master_matrix>
template <typename F>
inline void Chain_matrix<Master_matrix>::_add_to(Column& target, F&& addition)
{
  auto pivot = target.get_pivot();
  addition();

  if (pivot != target.get_pivot()) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      std::swap(pivotToColumnIndex_.at(pivot), pivotToColumnIndex_.at(target.get_pivot()));
    } else {
      std::swap(pivotToColumnIndex_[pivot], pivotToColumnIndex_[target.get_pivot()]);
    }
  }
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_remove_last(Index lastIndex)
{
  static_assert(Master_matrix::Option_list::has_removable_columns,
                "'_remove_last' is not implemented for the chosen options.");
  static_assert(Master_matrix::Option_list::has_map_column_container || !Master_matrix::Option_list::has_vine_update,
                "'_remove_last' is not implemented for the chosen options.");

  ID_index pivot;

  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    auto itToErase = matrix_.find(lastIndex);
    Column& colToErase = itToErase->second;
    pivot = colToErase.get_pivot();

    if constexpr (Master_matrix::Option_list::has_matrix_maximal_dimension_access) {
      Dim_opt::update_down(colToErase.get_dimension());
    }

    if (colToErase.is_paired()) matrix_.at(colToErase.get_paired_chain_index()).unassign_paired_chain();
    pivotToColumnIndex_.erase(pivot);
    matrix_.erase(itToErase);
  } else {
    GUDHI_CHECK(lastIndex == nextIndex_ - 1 && nextIndex_ == matrix_.size(),
                std::logic_error("Chain_matrix::_remove_last - Indexation problem."));

    Column& colToErase = matrix_[lastIndex];
    pivot = colToErase.get_pivot();

    if constexpr (Master_matrix::Option_list::has_matrix_maximal_dimension_access) {
      Dim_opt::update_down(colToErase.get_dimension());
    }

    if (colToErase.is_paired()) matrix_.at(colToErase.get_paired_chain_index()).unassign_paired_chain();
    pivotToColumnIndex_[pivot] = -1;
    matrix_.pop_back();
    // TODO: resize matrix_ when a lot is removed? Could be not the best strategy if user inserts a lot back afterwards.
  }

  if constexpr (!Master_matrix::Option_list::has_vine_update) {
    --nextIndex_;  // should not be updated when there are vine updates, as possibly lastIndex != nextIndex - 1
  }

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    auto it = _indexToBar().find(--_nextPosition());
    typename Barcode::iterator bar = it->second;

    if (bar->death == static_cast<Pos_index>(-1))
      _barcode().erase(bar);
    else
      bar->death = -1;

    _indexToBar().erase(it);
    if constexpr (Master_matrix::Option_list::has_vine_update) Swap_opt::CP::pivotToPosition_.erase(pivot);
  }

  if constexpr (Master_matrix::Option_list::has_row_access) {
    GUDHI_CHECK(
        RA_opt::get_row(pivot).size() == 0,
        std::invalid_argument(
            "Chain_matrix::_remove_last - Column asked to be removed does not corresponds to a maximal simplex."));
    if constexpr (Master_matrix::Option_list::has_removable_rows) {
      RA_opt::erase_empty_row(pivot);
    }
  }
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_update_barcode(Pos_index birth)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      auto& barIt = _indexToBar().at(birth);
      barIt->death = _nextPosition();
      _indexToBar().try_emplace(_nextPosition(), barIt);  // list so iterators are stable
    } else {
      _barcode()[_indexToBar()[birth]].death = _nextPosition();
      _indexToBar().push_back(_indexToBar()[birth]);
    }
    ++_nextPosition();
  }
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_add_bar(Dimension dim)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _barcode().emplace_back(_nextPosition(), -1, dim);
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      _indexToBar().try_emplace(_nextPosition(), --_barcode().end());
    } else {
      _indexToBar().push_back(_barcode().size() - 1);
    }
    ++_nextPosition();
  }
}

template <class Master_matrix>
template <class Container>
inline void Chain_matrix<Master_matrix>::_container_insert(const Container& column, Index pos, Dimension dim)
{
  ID_index pivot;
  if constexpr (Master_matrix::Option_list::is_z2) {
    pivot = *(column.rbegin());
  } else {
    pivot = column.rbegin()->first;
  }
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    pivotToColumnIndex_.try_emplace(pivot, pos);
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.try_emplace(pos, Column(pos, column, dim, RA_opt::rows_, colSettings_));
    } else {
      matrix_.try_emplace(pos, Column(column, dim, colSettings_));
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.emplace_back(pos, column, dim, RA_opt::rows_, colSettings_);
    } else {
      matrix_.emplace_back(column, dim, colSettings_);
    }
    pivotToColumnIndex_[pivot] = pos;
  }
}

template <class Master_matrix>
inline void Chain_matrix<Master_matrix>::_container_insert(const Column& column, [[maybe_unused]] Index pos)
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.try_emplace(pos, Column(column, column.get_column_index(), RA_opt::rows_, colSettings_));
    } else {
      matrix_.try_emplace(pos, Column(column, colSettings_));
    }
  } else {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      matrix_.emplace_back(column, column.get_column_index(), RA_opt::rows_, colSettings_);
    } else {
      matrix_.emplace_back(column, colSettings_);
    }
  }
}

template <class Master_matrix>
inline constexpr typename Chain_matrix<Master_matrix>::Barcode& Chain_matrix<Master_matrix>::_barcode()
{
  if constexpr (Master_matrix::Option_list::has_vine_update)
    return Swap_opt::template Chain_pairing<Master_matrix>::barcode_;
  else
    return Pair_opt::barcode_;
}

template <class Master_matrix>
inline constexpr typename Chain_matrix<Master_matrix>::Bar_dictionary&
Chain_matrix<Master_matrix>::_indexToBar()
{
  if constexpr (Master_matrix::Option_list::has_vine_update)
    return Swap_opt::template Chain_pairing<Master_matrix>::indexToBar_;
  else
    return Pair_opt::indexToBar_;
}

template <class Master_matrix>
inline constexpr typename Chain_matrix<Master_matrix>::Pos_index& Chain_matrix<Master_matrix>::_nextPosition()
{
  if constexpr (Master_matrix::Option_list::has_vine_update)
    return Swap_opt::template Chain_pairing<Master_matrix>::nextPosition_;
  else
    return Pair_opt::nextPosition_;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_CHAIN_MATRIX_H
