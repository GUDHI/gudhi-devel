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
 * @file Id_to_index_overlay.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Id_to_index_overlay class.
 */

#ifndef PM_ID_TO_POS_TRANSLATION_H
#define PM_ID_TO_POS_TRANSLATION_H

#include <cmath>
#include <vector>
#include <cassert>
#include <utility>    //std::swap, std::move & std::exchange
#include <algorithm>  //std::transform
#include <stdexcept>  //std::invalid_argument

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Id_to_index_overlay Id_to_index_overlay.h gudhi/Persistence_matrix/Id_to_index_overlay.h
 * @ingroup persistence_matrix
 *
 * @brief Overlay for @ref mp_matrices "non-basic matrices" replacing all input and output @ref MatIdx indices of
 * the original methods with @ref IDIdx indices.
 * 
 * @tparam Underlying_matrix %Matrix type taking the overlay.
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Underlying_matrix, class Master_matrix>
class Id_to_index_overlay 
{
 public:
  using Index = typename Master_matrix::Index;                          /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;                    /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index;                  /**< @ref PosIdx index type. */
  using Dimension = typename Master_matrix::Dimension;                  /**< Dimension value type. */
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element = typename Master_matrix::Element;                /**< Type of an field element. */
  using Boundary = typename Master_matrix::Boundary;                    /**< Type of an input column. */
  using Column = typename Master_matrix::Column;                        /**< Column type. */
  using Row = typename Master_matrix::Row;                              /**< Row type,
                                                                             only necessary with row access option. */
  using Bar = typename Master_matrix::Bar;                              /**< Bar type. */
  using Barcode = typename Master_matrix::Barcode;                      /**< Barcode type. */
  using Cycle = typename Master_matrix::Cycle;                          /**< Cycle type. */
  using Entry_constructor = typename Master_matrix::Entry_constructor;  /**< Factory of @ref Entry classes. */
  using Column_settings = typename Master_matrix::Column_settings;      /**< Structure giving access to the columns to
                                                                             necessary external classes. */

  /**
   * @brief Constructs an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Id_to_index_overlay(Column_settings* colSettings);
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
   * interest and not everything should be stored, then use the @ref insert_boundary method instead
   * (after creating the matrix with the @ref Id_to_index_overlay(unsigned int, Column_settings*)
   * constructor preferably).
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  template <class Boundary_range = Boundary>
  Id_to_index_overlay(const std::vector<Boundary_range>& orderedBoundaries, 
                      Column_settings* colSettings);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Id_to_index_overlay(unsigned int numberOfColumns, Column_settings* colSettings);
  /**
   * @brief Only available for @ref chainmatrix "chain matrices". Constructs an empty matrix and stores the given
   * comparators.
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
  Id_to_index_overlay(Column_settings* colSettings,
                      const BirthComparatorFunction& birthComparator, 
                      const DeathComparatorFunction& deathComparator);
  /**
   * @brief Only available for @ref chainmatrix "chain matrices". 
   * Constructs a new matrix from the given ranges of @ref Matrix::Entry_representative. Each range corresponds to a
   * column (the order of the ranges are preserved). The content of the ranges is assumed to be sorted by increasing
   * IDs. The IDs of the simplices are also assumed to be consecutive, ordered by filtration value, starting with 0. 
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
   * (after creating the matrix with the @ref Id_to_index_overlay(unsigned int, Column_settings*,
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
  template <typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_range>
  Id_to_index_overlay(const std::vector<Boundary_range>& orderedBoundaries, 
                      Column_settings* colSettings, 
                      const BirthComparatorFunction& birthComparator, 
                      const DeathComparatorFunction& deathComparator);
  /**
   * @brief Only available for @ref chainmatrix "chain matrices".
   * Constructs a new empty matrix and reserves space for the given number of columns.
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
  Id_to_index_overlay(unsigned int numberOfColumns, 
                      Column_settings* colSettings,
                      const BirthComparatorFunction& birthComparator, 
                      const DeathComparatorFunction& deathComparator);
  /**
   * @brief Copy constructor. If @p operators or @p entryConstructor is not a null pointer, its value is kept
   * instead of the one in the copied matrix.
   * 
   * @param matrixToCopy Matrix to copy.
   * @param colSettings Either a pointer to an existing setting structure for the columns or a null pointer.
   * The structure should contain all the necessary external classes specifically necessary for the choosen column type,
   * such as custom allocators. If null pointer, the pointer stored in @p matrixToCopy is used instead.
   */
  Id_to_index_overlay(const Id_to_index_overlay& matrixToCopy, 
                      Column_settings* colSettings = nullptr);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Id_to_index_overlay(Id_to_index_overlay&& other) noexcept;
  /**
   * @brief Destructor.
   */
  ~Id_to_index_overlay();

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
   * The content of the new column will vary depending on the underlying @ref mp_matrices "type of the matrix":
   * - If it is a boundary type matrix and only \f$ R \f$ is stored, the boundary is just copied. The column will only 
   *   be reduced later when the barcode is requested in order to apply some optimizations with the additional
   *   knowledge. Hence, the barcode will also not be updated.
   * - If it is a boundary type matrix and both \f$ R \f$ and \f$ U \f$ are stored, the new boundary is stored in its
   *   reduced form and the barcode, if active, is also updated.
   * - If it is a chain type matrix, the new column is of the form 
   *   `IDIdx + linear combination of older column IDIdxs`, where the combination is deduced while reducing the 
   *   given boundary. If the barcode is stored, it will also be updated.
   * 
   * @tparam Boundary_range Range of @ref Matrix::Entry_representative. Assumed to have a begin(), end() and size()
   * method.
   * @param boundary Boundary generating the new column. The content should be ordered by ID.
   * @param dim Dimension of the cell whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   */
  template <class Boundary_range = Boundary>
  void insert_boundary(const Boundary_range& boundary, Dimension dim = -1);
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
  void insert_boundary(ID_index cellIndex, const Boundary_range& boundary, Dimension dim = -1);
  /**
   * @brief Returns the column at the given @ref IDIdx index. 
   * For @ref boundarymatrix "RU matrices", the returned column is from \f$ R \f$.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param cellID @ref IDIdx index of the column to return.
   * @return Reference to the column.
   */
  Column& get_column(ID_index cellID);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access is true.
   * Returns the row at the given @ref rowindex "row index".
   * For @ref boundarymatrix "RU matrices", the returned row is from \f$ R \f$.
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   *
   * @warning The @ref Entry_column_index::get_column_index "get_column_index" method of the row entries returns the
   * original @ref PosIdx indices (before any swaps) for @ref boundarymatrix "boundary matrices" and
   * @ref MatIdx indices for @ref chainmatrix "chain matrices".
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return: @ref IDIdx for @ref chainmatrix "chain matrices"
   * or updated @ref IDIdx for @ref boundarymatrix "boundary matrices" if swaps occurred.
   * @return Reference to the row.
   */
  Row& get_row(ID_index rowIndex);
  /**
   * @brief The effect varies depending on the matrices and the options:
   * - @ref boundarymatrix "boundary matrix" with only \f$ R \f$ stored:
   *    - @ref PersistenceMatrixOptions::has_map_column_container and
   *      @ref PersistenceMatrixOptions::has_column_and_row_swaps are true:
   *      cleans up maps used for the lazy row swaps.
   *    - @ref PersistenceMatrixOptions::has_row_access and @ref PersistenceMatrixOptions::has_removable_rows are true:
   *      assumes that the row is empty and removes it. 
   *    - Otherwise, does nothing.
   * - @ref boundarymatrix "boundary matrix" with \f$ U \f$ stored: only \f$ R \f$ is affected by the above.
   *   If properly used, \f$ U \f$ will never have empty rows.
   * - @ref chainmatrix "chain matrix": only available if @ref PersistenceMatrixOptions::has_row_access and
   *   @ref PersistenceMatrixOptions::has_removable_rows are true.
   *   Assumes that the row is empty and removes it. 
   *
   * @warning The removed rows are always assumed to be empty. If it is not the case, the deleted row entries are not
   * removed from their columns. And in the case of intrusive rows, this will generate a segmentation fault when 
   * the column entries are destroyed later. The row access is just meant as a "read only" access to the rows and the
   * @ref erase_empty_row method just as a way to specify that a row is empty and can therefore be removed from
   * dictionaries. This allows to avoid testing the emptiness of a row at each column entry removal, what can be
   * quite frequent. 
   * 
   * @param rowIndex @ref rowindex "Row index" of the empty row to remove.
   */
  void erase_empty_row(ID_index rowIndex);
  /**
   * @brief Only available for RU and @ref chainmatrix "chain matrices" and if
   * @ref PersistenceMatrixOptions::has_removable_columns and @ref PersistenceMatrixOptions::has_vine_update are true.
   * For @ref chainmatrix "chain matrices", @ref PersistenceMatrixOptions::has_map_column_container and
   * @ref PersistenceMatrixOptions::has_column_pairings also need to be true.
   * Assumes that the cell is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., RU is still an upper triangular decomposition of the @ref boundarymatrix "boundary matrix" and chain is
   * still a compatible bases of the chain complex in the sense of @cite zigzag).
   * The maximality of the cell is not verified.
   * Also updates the barcode if it was computed.
   *
   * For @ref chainmatrix "chain matrices", using the other version of the method could perform better depending on
   * how the data is maintained on the side of the user. Then, @ref PersistenceMatrixOptions::has_column_pairings also
   * do not need to be true.
   *
   * See also @ref remove_last.
   * 
   * @param cellID @ref IDIdx index of the cell to remove.
   */
  void remove_maximal_cell(ID_index cellID);
  /**
   * @brief Only available for @ref chainmatrix "chain matrices" and if
   * @ref PersistenceMatrixOptions::has_removable_columns, @ref PersistenceMatrixOptions::has_vine_update
   * and @ref PersistenceMatrixOptions::has_map_column_container are true.
   * Assumes that the cell is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., it is still a compatible bases of the chain complex in the sense of @cite zigzag).
   * The maximality of the cell is not verified.
   * Also updates the barcode if it was computed.
   *
   * To maintain the compatibility, vine swaps are done to move the cell up to the end of the filtration. Once at 
   * the end, the removal is trivial. But for @ref chainmatrix "chain matrices", swaps do not actually swap the position
   * of the column every time, so the cells appearing after @p cellIndex in the filtration have to be searched first
   * within the matrix. If the user has an easy access to the @ref IDIdx of the cells in the order of filtration,
   * passing them by argument with @p columnsToSwap allows to skip a linear search process. Typically, if the user knows
   * that the cell he wants to remove is already the last cell of the filtration, calling
   * @ref remove_maximal_cell(ID_index, const std::vector<ID_index>&) "remove_maximal_cell(cellID, {})"
   * will be faster than @ref remove_last().
   *
   * See also @ref remove_last.
   * 
   * @param cellID @ref IDIdx index of the cell to remove.
   * @param columnsToSwap Vector of @ref IDIdx indices of the cells coming after @p cellID in the filtration.
   */
  void remove_maximal_cell(ID_index cellID, const std::vector<ID_index>& columnsToSwap);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns is true. Additionally, if the
   * matrix is a @ref chainmatrix "chain matrix", either @ref PersistenceMatrixOptions::has_map_column_container has to
   * be true or @ref PersistenceMatrixOptions::has_vine_update has to be false.
   * Removes the last cell in the filtration from the matrix and updates the barcode if it is stored.
   * 
   * See also @ref remove_maximal_cell.
   *
   * For @ref chainmatrix "chain matrices", if @ref PersistenceMatrixOptions::has_vine_update is true, the last cell
   * does not have to be at the end of the matrix and therefore has to be searched first. In this case, if the user
   * already knows the @ref IDIdx of the last cell, calling
   * @ref remove_maximal_cell(ID_index, const std::vector<ID_index>&) "remove_maximal_cell(cellID, {})"
   * instead allows to skip the search.
   */
  void remove_last();

  /**
   * @brief Returns the maximal dimension of a cell stored in the matrix. Only available 
   * if @ref PersistenceMatrixOptions::has_matrix_maximal_dimension_access is true.
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
   * @brief Returns the dimension of the given cell. Only available for @ref mp_matrices "non-basic matrices".
   * 
   * @param cellID @ref IDIdx index of the cell.
   * @return Dimension of the cell.
   */
  Dimension get_column_dimension(ID_index cellID) const;

  /**
   * @brief Adds column corresponding to @p sourceCellID onto the column corresponding to @p targetCellID.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of the matrix.
   * For example, a right-to-left addition could corrupt the computation of the barcode if done blindly.
   * So should be used with care.
   * 
   * @param sourceCellID @ref IDIdx index of the source column.
   * @param targetCellID @ref IDIdx index of the target column.
   */
  void add_to(ID_index sourceCellID, ID_index targetCellID);
  /**
   * @brief Multiplies the target column with the coefficient and then adds the source column to it.
   * That is: `targetColumn = (targetColumn * coefficient) + sourceColumn`.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of the matrix.
   * For example, a right-to-left addition could corrupt the computation of the barcode if done blindly.
   * So should be used with care.
   * 
   * @param sourceCellID @ref IDIdx index of the source column.
   * @param coefficient Value to multiply.
   * @param targetCellID @ref IDIdx index of the target column.
   */
  void multiply_target_and_add_to(ID_index sourceCellID, const Field_element& coefficient, ID_index targetCellID);
  /**
   * @brief Multiplies the source column with the coefficient before adding it to the target column.
   * That is: `targetColumn += (coefficient * sourceColumn)`. The source column will **not** be modified.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of the matrix.
   * For example, a right-to-left addition could corrupt the computation of the barcode if done blindly.
   * So should be used with care.
   * 
   * @param coefficient Value to multiply.
   * @param sourceCellID @ref IDIdx index of the source column.
   * @param targetCellID @ref IDIdx index of the target column.
   */
  void multiply_source_and_add_to(const Field_element& coefficient, ID_index sourceCellID, ID_index targetCellID);

  /**
   * @brief Zeroes the entry at the given coordinates. Not available for @ref chainmatrix "chain matrices".
   * In general, should be used with care to not destroy the validity 
   * of the persistence related properties of the matrix.
   *
   * For @ref boundarymatrix "RU matrices", zeros only the entry in \f$ R \f$.
   * 
   * @param cellID @ref IDIdx index of the cell corresponding to the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   */
  void zero_entry(ID_index cellID, ID_index rowIndex);
  /**
   * @brief Zeroes the column at the given index. Not available for @ref chainmatrix "chain matrices".
   * In general, should be used with care to not destroy the validity 
   * of the persistence related properties of the matrix.
   *
   * For @ref boundarymatrix "RU matrices", zeros only the column in \f$ R \f$.
   * 
   * @param cellID @ref IDIdx index of the cell corresponding to the column.
   */
  void zero_column(ID_index cellID);
  /**
   * @brief Indicates if the entry at given coordinates has value zero.
   *
   * For @ref boundarymatrix "RU matrices", looks into \f$ R \f$.
   * 
   * @param cellID @ref IDIdx index of the cell corresponding to the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @return true If the entry has value zero.
   * @return false Otherwise.
   */
  bool is_zero_entry(ID_index cellID, ID_index rowIndex) const;
  /**
   * @brief Indicates if the column at given index has value zero.
   *
   * For @ref boundarymatrix "RU matrices", looks into \f$ R \f$.
   *
   * Note that for @ref chainmatrix "chain matrices", this method should always return false, as a valid
   * @ref chainmatrix "chain matrix" never has empty columns.
   * 
   * @param cellID @ref IDIdx index of the cell corresponding to the column.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(ID_index cellID);

  /**
   * @brief Returns the @ref IDIdx index of the column which has the given @ref rowindex "row index" as pivot.
   * Assumes that the pivot exists. For @ref boundarymatrix "RU matrices", the column is returned from \f$ R \f$.
   *
   * Recall that the row indices for @ref chainmatrix "chain matrices" correspond to the @ref IDIdx indices and that
   * the row indices for a @ref boundarymatrix "RU matrix" correspond to the updated @ref IDIdx indices which got
   * potentially swapped by a vine swap.
   * 
   * @param cellIndex @ref rowindex "Row index" of the pivot.
   * @return @ref IDIdx index of the column with the given pivot.
   */
  ID_index get_column_with_pivot(ID_index cellIndex) const;
  /**
   * @brief Returns the @ref rowindex "row index" of the pivot of the given column.
   * 
   * @param cellID @ref IDIdx index of the cell corresponding to the column.
   * @return The @ref rowindex "row index" of the pivot.
   */
  ID_index get_pivot(ID_index cellID);

  /**
   * @brief Resets the matrix to an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  void reset(Column_settings* colSettings) {
    matrix_.reset(colSettings);
    nextIndex_ = 0;
  }

  // void set_operators(Field_operators* operators) { matrix_.set_operators(operators); }

  /**
   * @brief Assign operator.
   */
  Id_to_index_overlay& operator=(const Id_to_index_overlay& other);
  /**
   * @brief Swap operator.
   */
  friend void swap(Id_to_index_overlay& matrix1, Id_to_index_overlay& matrix2) {
    swap(matrix1.matrix_, matrix2.matrix_);
    if (Master_matrix::Option_list::is_of_boundary_type) std::swap(matrix1.idToIndex_, matrix2.idToIndex_);
    std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
  }

  void print();  // for debug

  // access to optional methods

  /**
   * @brief Returns the current barcode of the matrix.
   * Available only if @ref PersistenceMatrixOptions::has_column_pairings is true.
   *
   * Recall that we assume that the boundaries were inserted in the order of filtration for the barcode to be valid.
   *
   * @warning For simple @ref boundarymatrix "boundary matrices" (only storing \f$ R \f$), we assume that
   * @ref get_current_barcode is only called once, when the matrix is completed.
   * 
   * @return A reference to the barcode. The barcode is a vector of @ref Matrix::Bar. A bar stores three informations:
   * the @ref PosIdx birth index, the @ref PosIdx death index and the dimension of the bar.
   */
  const Barcode& get_current_barcode();
  
  /**
   * @brief Only available for simple @ref boundarymatrix "boundary matrices" (only storing \f$ R \f$) and if
   * @ref PersistenceMatrixOptions::has_column_and_row_swaps is true.
   * Swaps the two given columns. Note that it really just swaps two columns and do not updates
   * anything else, nor performs additions to maintain some properties on the matrix.
   * 
   * @param cellID1 First column @ref IDIdx index to swap.
   * @param cellID2 Second column @ref IDIdx index to swap.
   */
  void swap_columns(ID_index cellID1, ID_index cellID2);
  /**
   * @brief Only available for simple @ref boundarymatrix "boundary matrices" (only storing R)
   * and if @ref PersistenceMatrixOptions::has_column_and_row_swaps is true.
   * Swaps the two given rows. Note that it really just swaps two rows and do not updates
   * anything else, nor performs additions to maintain some properties on the matrix.
   * 
   * @param rowIndex1 First @ref rowindex "row index" to swap.
   * @param rowIndex2 Second @ref rowindex "row index" to swap.
   */
  void swap_rows(Index rowIndex1, Index rowIndex2);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true.
   * Does the same than @ref vine_swap, but assumes that the swap is non trivial and
   * therefore skips a part of the case study.
   * 
   * @param cellID1 @ref IDIdx index of the first cell.
   * @param cellID2 @ref IDIdx index of the second cell. It is assumed that the @ref PosIdx of both only differs by one.
   * @return Let \f$ pos1 \f$ be the @ref PosIdx index of @p columnIndex1 and \f$ pos2 \f$ be the @ref PosIdx index of
   * @p columnIndex2. The method returns the @ref MatIdx of the column which has now, after the swap, the @ref PosIdx
   * \f$ max(pos1, pos2) \f$.
   */
  ID_index vine_swap_with_z_eq_1_case(ID_index cellID1, ID_index cellID2);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true.
   * Does a vine swap between two cells which are consecutive in the filtration. Roughly, if \f$ F \f$ is
   * the current filtration represented by the matrix, the method modifies the matrix such that the new state
   * corresponds to a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but with the two given cells
   * at swapped positions. Of course, the two cells should not have a face/coface relation which each other ;
   * \f$ F' \f$ has to be a valid filtration.
   * See @cite vineyards for more information about vine and vineyards.
   * 
   * @param cellID1 @ref IDIdx index of the first cell.
   * @param cellID2 @ref IDIdx index of the second cell. It is assumed that the @ref PosIdx of both only differs by one.
   * @return Let \f$ pos1 \f$ be the @ref PosIdx index of @p columnIndex1 and \f$ pos2 \f$ be the @ref PosIdx index of
   * @p columnIndex2. The method returns the @ref MatIdx of the column which has now, after the swap, the @ref PosIdx
   * \f$ max(pos1, pos2) \f$.
   */
  ID_index vine_swap(ID_index cellID1, ID_index cellID2);
  
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true. Pre-computes
   * the representative cycles of the current state of the filtration represented by the matrix.
   * It does not need to be called before `get_representative_cycles` is called for the first time, but needs to be
   * called before calling `get_representative_cycles` again if the matrix was modified in between. Otherwise the
   * old cycles will be returned.
   */
  void update_representative_cycles();
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true.
   * Returns all representative cycles of the current filtration.
   * 
   * @return A const reference to the vector of representative cycles.
   */
  const std::vector<Cycle>& get_representative_cycles();
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true.
   * Returns the cycle representing the given bar.
   * 
   * @param bar A bar from the current barcode.
   * @return A const reference to the cycle representing @p bar.
   */
  const Cycle& get_representative_cycle(const Bar& bar);

 private:
  using Dictionary = typename Master_matrix::template Dictionary<Index>;

  Underlying_matrix matrix_;  /**< Interfaced matrix. */
  Dictionary* idToIndex_;     /**< Map from @ref IDIdx index to @ref MatIdx index. */
  Index nextIndex_;           /**< Next unused index. */

  void _initialize_map(unsigned int size);
  Index _id_to_index(ID_index id) const;
  Index& _id_to_index(ID_index id);
};

template <class Underlying_matrix, class Master_matrix>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>::Id_to_index_overlay(Column_settings* colSettings)
    : matrix_(colSettings), idToIndex_(nullptr), nextIndex_(0) 
{
  _initialize_map(0);
}

template <class Underlying_matrix, class Master_matrix>
template <class Boundary_range>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>::Id_to_index_overlay(
    const std::vector<Boundary_range>& orderedBoundaries, Column_settings* colSettings)
    : matrix_(orderedBoundaries, colSettings), idToIndex_(nullptr), nextIndex_(orderedBoundaries.size()) 
{
  _initialize_map(orderedBoundaries.size());
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    for (unsigned int i = 0; i < orderedBoundaries.size(); i++) {
      _id_to_index(i) = i;
    }
  }
}

template <class Underlying_matrix, class Master_matrix>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>::Id_to_index_overlay(unsigned int numberOfColumns,
                                                                                 Column_settings* colSettings)
    : matrix_(numberOfColumns, colSettings), idToIndex_(nullptr), nextIndex_(0) 
{
  _initialize_map(numberOfColumns);
}

template <class Underlying_matrix, class Master_matrix>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>::Id_to_index_overlay(
    Column_settings* colSettings, 
    const BirthComparatorFunction& birthComparator, 
    const DeathComparatorFunction& deathComparator)
    : matrix_(colSettings, birthComparator, deathComparator), idToIndex_(nullptr), nextIndex_(0) 
{
  _initialize_map(0);
}

template <class Underlying_matrix, class Master_matrix>
template <typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_range>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>::Id_to_index_overlay(
    const std::vector<Boundary_range>& orderedBoundaries, 
    Column_settings* colSettings,
    const BirthComparatorFunction& birthComparator, 
    const DeathComparatorFunction& deathComparator)
    : matrix_(orderedBoundaries, colSettings, birthComparator, deathComparator),
      idToIndex_(nullptr),
      nextIndex_(orderedBoundaries.size()) 
{
  _initialize_map(orderedBoundaries.size());
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    for (unsigned int i = 0; i < orderedBoundaries.size(); i++) {
      _id_to_index(i) = i;
    }
  }
}

template <class Underlying_matrix, class Master_matrix>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>::Id_to_index_overlay(
    unsigned int numberOfColumns, 
    Column_settings* colSettings,
    const BirthComparatorFunction& birthComparator, 
    const DeathComparatorFunction& deathComparator)
    : matrix_(numberOfColumns, colSettings, birthComparator, deathComparator),
      idToIndex_(nullptr),
      nextIndex_(0) 
{
  _initialize_map(numberOfColumns);
}

template <class Underlying_matrix, class Master_matrix>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>::Id_to_index_overlay(
    const Id_to_index_overlay& matrixToCopy, Column_settings* colSettings)
    : matrix_(matrixToCopy.matrix_, colSettings),
      idToIndex_(nullptr),
      nextIndex_(matrixToCopy.nextIndex_) 
{
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    idToIndex_ = new Dictionary(*matrixToCopy.idToIndex_);
  } else {
    idToIndex_ = &matrix_.pivotToColumnIndex_;
  }
}

template <class Underlying_matrix, class Master_matrix>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>::Id_to_index_overlay(Id_to_index_overlay&& other) noexcept
    : matrix_(std::move(other.matrix_)),
      idToIndex_(std::exchange(other.idToIndex_, nullptr)),
      nextIndex_(std::exchange(other.nextIndex_, 0)) 
{}

template <class Underlying_matrix, class Master_matrix>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>::~Id_to_index_overlay() 
{
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    if (idToIndex_ != nullptr) delete idToIndex_;
  }
}

template <class Underlying_matrix, class Master_matrix>
template <class Boundary_range>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::insert_boundary(const Boundary_range& boundary,
                                                                                  Dimension dim) 
{
  matrix_.insert_boundary(boundary, dim);
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      idToIndex_->emplace(nextIndex_, nextIndex_);
    } else {
      if (idToIndex_->size() == nextIndex_) {
        idToIndex_->push_back(nextIndex_);
      } else {
        _id_to_index(nextIndex_) = nextIndex_;
      }
    }
    ++nextIndex_;
  }
}

template <class Underlying_matrix, class Master_matrix>
template <class Boundary_range>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::insert_boundary(ID_index cellIndex,
                                                                                  const Boundary_range& boundary,
                                                                                  Dimension dim) 
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    GUDHI_CHECK(idToIndex_->find(cellIndex) == idToIndex_->end(),
                std::invalid_argument("Id_to_index_overlay::insert_boundary - Index for simplex already chosen!"));
  } else {
    GUDHI_CHECK((idToIndex_->size() <= cellIndex || _id_to_index(cellIndex) == static_cast<Index>(-1)),
                std::invalid_argument("Id_to_index_overlay::insert_boundary - Index for simplex already chosen!"));
  }
  matrix_.insert_boundary(cellIndex, boundary, dim);
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      idToIndex_->emplace(cellIndex, nextIndex_);
    } else {
      if (idToIndex_->size() <= cellIndex) {
        idToIndex_->resize(cellIndex + 1, -1);
      }
      _id_to_index(cellIndex) = nextIndex_;
    }
    ++nextIndex_;
  }
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Column&
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_column(ID_index cellID) 
{
  return matrix_.get_column(_id_to_index(cellID));
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Row&
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_row(ID_index rowIndex) 
{
  return matrix_.get_row(rowIndex);
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::erase_empty_row(ID_index rowIndex) 
{
  return matrix_.erase_empty_row(rowIndex);
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::remove_maximal_cell(ID_index cellID) 
{
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    std::vector<ID_index> indexToID(nextIndex_);
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      for (auto& p : *idToIndex_) {
        indexToID[p.second] = p.first;
      }
    } else {
      for (ID_index i = 0; i < idToIndex_->size(); ++i) {
        if (_id_to_index(i) != static_cast<Index>(-1)) indexToID[_id_to_index(i)] = i;
      }
    }
    --nextIndex_;
    for (Index curr = _id_to_index(cellID); curr < nextIndex_; ++curr) {
      matrix_.vine_swap(curr);
      std::swap(idToIndex_->at(indexToID[curr]), idToIndex_->at(indexToID[curr + 1]));
    }
    matrix_.remove_last();
    GUDHI_CHECK(_id_to_index(cellID) == nextIndex_,
                std::logic_error("Id_to_index_overlay::remove_maximal_cell - Indexation problem."));

    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      idToIndex_->erase(cellID);
    } else {
      _id_to_index(cellID) = -1;
    }
  } else {
    matrix_.remove_maximal_cell(cellID);
  }
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::remove_maximal_cell(
    ID_index cellID, const std::vector<ID_index>& columnsToSwap) 
{
  static_assert(!Master_matrix::Option_list::is_of_boundary_type,
                "'remove_maximal_cell(ID_index,const std::vector<Index>&)' is not available for the chosen options.");
  std::vector<Index> translatedIndices;
  std::transform(columnsToSwap.cbegin(), columnsToSwap.cend(), std::back_inserter(translatedIndices),
                 [&](ID_index id) { return _id_to_index(id); });
  matrix_.remove_maximal_cell(cellID, translatedIndices);
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::remove_last() 
{
  if (idToIndex_->empty()) return;  //empty matrix

  matrix_.remove_last();

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    --nextIndex_;
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      auto it = idToIndex_->begin();
      while (it->second != nextIndex_) ++it;   //should never reach idToIndex_->end()
      idToIndex_->erase(it);
    } else {
      Index id = idToIndex_->size() - 1;
      while (_id_to_index(id) == static_cast<Index>(-1)) --id;  // should always stop before reaching -1
      GUDHI_CHECK(_id_to_index(id) == nextIndex_,
                  std::logic_error("Id_to_index_overlay::remove_last - Indexation problem."));
      _id_to_index(id) = -1;
    }
  }
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Dimension
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_max_dimension() const 
{
  return matrix_.get_max_dimension();
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Index
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_number_of_columns() const 
{
  return matrix_.get_number_of_columns();
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Dimension
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_column_dimension(ID_index cellID) const 
{
  return matrix_.get_column_dimension(_id_to_index(cellID));
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::add_to(ID_index sourceCellID, ID_index targetCellID) 
{
  return matrix_.add_to(_id_to_index(sourceCellID), _id_to_index(targetCellID));
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::multiply_target_and_add_to(
    ID_index sourceCellID, const Field_element& coefficient, ID_index targetCellID) 
{
  return matrix_.multiply_target_and_add_to(_id_to_index(sourceCellID), coefficient, _id_to_index(targetCellID));
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::multiply_source_and_add_to(
    const Field_element& coefficient, ID_index sourceCellID, ID_index targetCellID) 
{
  return matrix_.multiply_source_and_add_to(coefficient, _id_to_index(sourceCellID), _id_to_index(targetCellID));
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::zero_entry(ID_index cellID, ID_index rowIndex) 
{
  return matrix_.zero_entry(_id_to_index(cellID), rowIndex);
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::zero_column(ID_index cellID) 
{
  return matrix_.zero_column(_id_to_index(cellID));
}

template <class Underlying_matrix, class Master_matrix>
inline bool Id_to_index_overlay<Underlying_matrix, Master_matrix>::is_zero_entry(ID_index cellID,
                                                                               ID_index rowIndex) const 
{
  return matrix_.is_zero_entry(_id_to_index(cellID), rowIndex);
}

template <class Underlying_matrix, class Master_matrix>
inline bool Id_to_index_overlay<Underlying_matrix, Master_matrix>::is_zero_column(ID_index cellID) 
{
  return matrix_.is_zero_column(_id_to_index(cellID));
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::ID_index
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_column_with_pivot(ID_index simplexIndex) const 
{
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    Index pos = matrix_.get_column_with_pivot(simplexIndex);
    ID_index i = 0;
    while (_id_to_index(i) != pos) ++i;
    return i;
  } else {
    return simplexIndex;
  }
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::ID_index
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_pivot(ID_index cellID) 
{
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    return matrix_.get_pivot(_id_to_index(cellID));
  } else {
    return cellID;
  }
}

template <class Underlying_matrix, class Master_matrix>
inline Id_to_index_overlay<Underlying_matrix, Master_matrix>&
Id_to_index_overlay<Underlying_matrix, Master_matrix>::operator=(const Id_to_index_overlay& other) 
{
  matrix_ = other.matrix_;
  if (Master_matrix::Option_list::is_of_boundary_type)
    idToIndex_ = other.idToIndex_;
  else
    idToIndex_ = &matrix_.pivotToColumnIndex_;
  nextIndex_ = other.nextIndex_;

  return *this;
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::print() 
{
  return matrix_.print();
}

template <class Underlying_matrix, class Master_matrix>
inline const typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Barcode&
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_current_barcode() 
{
  return matrix_.get_current_barcode();
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::update_representative_cycles() 
{
  matrix_.update_representative_cycles();
}

template <class Underlying_matrix, class Master_matrix>
inline const std::vector<typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Cycle>&
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_representative_cycles() 
{
  return matrix_.get_representative_cycles();
}

template <class Underlying_matrix, class Master_matrix>
inline const typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Cycle&
Id_to_index_overlay<Underlying_matrix, Master_matrix>::get_representative_cycle(const Bar& bar) 
{
  return matrix_.get_representative_cycle(bar);
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::swap_columns(ID_index cellID1, ID_index cellID2) 
{
  matrix_.swap_columns(_id_to_index(cellID1), _id_to_index(cellID2));
  std::swap(idToIndex_->at(cellID1), idToIndex_->at(cellID2));
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::swap_rows(Index rowIndex1, Index rowIndex2) 
{
  matrix_.swap_rows(rowIndex1, rowIndex2);
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::ID_index
Id_to_index_overlay<Underlying_matrix, Master_matrix>::vine_swap_with_z_eq_1_case(ID_index cellID1, ID_index cellID2) 
{
  Index first = _id_to_index(cellID1);
  Index second = _id_to_index(cellID2);
  if (first > second) std::swap(first, second);

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    GUDHI_CHECK(second - first == 1,
                std::invalid_argument(
                    "Id_to_index_overlay::vine_swap_with_z_eq_1_case - The columns to swap are not contiguous."));

    bool change = matrix_.vine_swap_with_z_eq_1_case(first);

    std::swap(idToIndex_->at(cellID1), idToIndex_->at(cellID2));

    if (change) {
      return cellID1;
    }
    return cellID2;
  } else {
    return matrix_.vine_swap_with_z_eq_1_case(first, second);
  }
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::ID_index
Id_to_index_overlay<Underlying_matrix, Master_matrix>::vine_swap(ID_index cellID1, ID_index cellID2) 
{
  Index first = _id_to_index(cellID1);
  Index second = _id_to_index(cellID2);
  if (first > second) std::swap(first, second);

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    GUDHI_CHECK(second - first == 1,
                std::invalid_argument("Id_to_index_overlay::vine_swap - The columns to swap are not contiguous."));

    bool change = matrix_.vine_swap(first);

    std::swap(idToIndex_->at(cellID1), idToIndex_->at(cellID2));

    if (change) {
      return cellID1;
    }
    return cellID2;
  } else {
    return matrix_.vine_swap(first, second);
  }
}

template <class Underlying_matrix, class Master_matrix>
inline void Id_to_index_overlay<Underlying_matrix, Master_matrix>::_initialize_map([[maybe_unused]] unsigned int size) 
{
  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      idToIndex_ = new Dictionary(size);
    } else {
      idToIndex_ = new Dictionary(size, -1);
    }
  } else {
    idToIndex_ = &matrix_.pivotToColumnIndex_;
  }
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Index
Id_to_index_overlay<Underlying_matrix, Master_matrix>::_id_to_index(ID_index id) const 
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    return idToIndex_->at(id);
  } else {
    return idToIndex_->operator[](id);
  }
}

template <class Underlying_matrix, class Master_matrix>
inline typename Id_to_index_overlay<Underlying_matrix, Master_matrix>::Index&
Id_to_index_overlay<Underlying_matrix, Master_matrix>::_id_to_index(ID_index id)
{
  return idToIndex_->operator[](id);  //for maps, the entry is created if not existing as needed in the constructors
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_ID_TO_POS_TRANSLATION_H
