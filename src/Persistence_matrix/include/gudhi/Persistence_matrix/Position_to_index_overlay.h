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
 * @file Position_to_index_overlay.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Position_to_index_overlay class.
 */

#ifndef PM_POS_TO_ID_TRANSLATION_H
#define PM_POS_TO_ID_TRANSLATION_H

#include <vector>
#include <utility>    //std::swap, std::move & std::exchange
#include <algorithm>  //std::transform

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Position_to_index_overlay Position_to_index_overlay.h gudhi/Persistence_matrix/Position_to_index_overlay.h
 * @ingroup persistence_matrix
 *
 * @brief Overlay for @ref chainmatrix "chain matrices" replacing all input and output @ref MatIdx indices of the
 * original methods with @ref PosIdx indices. The overlay is useless for @ref boundarymatrix "boundary matrices"
 * as @ref MatIdx == @ref PosIdx for them.
 * 
 * @tparam %Underlying_matrix Matrix type taking the overlay.
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Underlying_matrix, class Master_matrix>
class Position_to_index_overlay 
{
 public:
  using Index = typename Master_matrix::Index;                              /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;                        /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index;                      /**< @ref PosIdx index type. */
  using Dimension = typename Master_matrix::Dimension;                      /**< Dimension value type. */
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element = typename Master_matrix::Element;                    /**< Type of an field element. */
  using Boundary = typename Master_matrix::Boundary;                        /**< Type of an input column. */
  using Column = typename Master_matrix::Column;                            /**< Column type. */
  using Row = typename Master_matrix::Row;                                  /**< Row type, only
                                                                                 necessary with row access option. */
  using Bar = typename Master_matrix::Bar;                                  /**< Bar type. */
  using Barcode = typename Master_matrix::Barcode;                          /**< Barcode type. */
  using Cycle = typename Master_matrix::Cycle;                              /**< Cycle type. */
  using Entry_representative = typename Master_matrix::Entry_representative;  /**< %Entry content representative. */
  using Entry_constructor = typename Master_matrix::Entry_constructor;        /**< Factory of @ref Entry classes. */
  using Column_settings = typename Master_matrix::Column_settings;          /**< Structure giving access to the columns
                                                                                 to necessary external classes. */

  /**
   * @brief Constructs an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Position_to_index_overlay(Column_settings* colSettings);
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
   * matrix with the @ref Position_to_index_overlay(unsigned int, Column_settings*)
   * constructor preferably).
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  template <class Boundary_range = Boundary>
  Position_to_index_overlay(const std::vector<Boundary_range>& orderedBoundaries, 
                            Column_settings* colSettings);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Position_to_index_overlay(unsigned int numberOfColumns, 
                            Column_settings* colSettings);
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
  Position_to_index_overlay(Column_settings* colSettings,
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
   * (after creating the matrix with the @ref Position_to_index_overlay(unsigned int, Column_settings*,
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
  Position_to_index_overlay(const std::vector<Boundary_range>& orderedBoundaries, 
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
  Position_to_index_overlay(unsigned int numberOfColumns, 
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
  Position_to_index_overlay(const Position_to_index_overlay& matrixToCopy, 
                            Column_settings* colSettings = nullptr);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Position_to_index_overlay(Position_to_index_overlay&& other) noexcept;

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
   * When inserted, the given boundary is reduced and from the reduction process, the column is deduced in the form of:
   * `IDIdx + linear combination of older column IDIdxs`. If the barcode is stored, it will be updated.
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
   * @p cellID values of precedent calls of the method for the corresponding cells and should be ordered in 
   * increasing order.
   * @param dim Dimension of the cell whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   */
  template <class Boundary_range = Boundary>
  void insert_boundary(ID_index cellIndex, const Boundary_range& boundary, Dimension dim = -1);
  /**
   * @brief Returns the column at the given @ref PosIdx index.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param position @ref PosIdx index of the column to return.
   * @return Reference to the column.
   */
  Column& get_column(Pos_index position);
  /**
   * @brief Returns the column at the given @ref PosIdx index.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param position @ref PosIdx index of the column to return.
   * @return Const reference to the column.
   */
  const Column& get_column(Pos_index position) const;
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access is true.
   * Returns the row at the given @ref rowindex "row index".
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return.
   * @return Reference to the row.
   */
  Row& get_row(ID_index rowIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access is true.
   * Returns the row at the given @ref rowindex "row index".
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return.
   * @return Const reference to the row.
   */
  const Row& get_row(ID_index rowIndex) const;
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access and
   * @ref PersistenceMatrixOptions::has_removable_rows are true.
   * Assumes that the row is empty and removes it. 
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
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns,
   * @ref PersistenceMatrixOptions::has_vine_update and @ref PersistenceMatrixOptions::has_map_column_container
   * are true.
   * Assumes that the cell is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., the matrix is still a compatible bases of the chain complex in the sense of @cite zigzag).
   * The maximality of the cell is not verified.
   * Also updates the barcode if it was computed.
   *
   * See also @ref remove_last.
   * 
   * @param position @ref PosIdx index of the cell to remove.
   */
  void remove_maximal_cell(Pos_index position);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns is true and,
   * if @ref PersistenceMatrixOptions::has_map_column_container is true or
   * @ref PersistenceMatrixOptions::has_vine_update is false.
   * Removes the last cell in the filtration from the matrix and updates the barcode if it is stored.
   *
   * See also @ref remove_maximal_cell.
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
   * @brief Returns the dimension of the given cell.
   * 
   * @param position @ref PosIdx index of the cell.
   * @return Dimension of the cell.
   */
  Dimension get_column_dimension(Pos_index position) const;

  /**
   * @brief Adds column corresponding to @p sourcePosition onto the column corresponding to @p targetPosition.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of the matrix.
   * For example, a right-to-left addition could corrupt the computation of the barcode if done blindly.
   * So should be used with care.
   * 
   * @param sourcePosition @ref PosIdx index of the source column.
   * @param targetPosition @ref PosIdx index of the target column.
   */
  void add_to(Pos_index sourcePosition, Pos_index targetPosition);
  /**
   * @brief Multiplies the target column with the coefficient and then adds the source column to it.
   * That is: `targetColumn = (targetColumn * coefficient) + sourceColumn`.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of the matrix.
   * For example, a right-to-left addition could corrupt the computation of the barcode if done blindly.
   * So should be used with care.
   * 
   * @param sourcePosition @ref PosIdx index of the source column.
   * @param coefficient Value to multiply.
   * @param targetPosition @ref PosIdx index of the target column.
   */
  void multiply_target_and_add_to(Pos_index sourcePosition, 
                                  const Field_element& coefficient,
                                  Pos_index targetPosition);
  /**
   * @brief Multiplies the source column with the coefficient before adding it to the target column.
   * That is: `targetColumn += (coefficient * sourceColumn)`. The source column will **not** be modified.
   *
   * @warning They will be no verification to ensure that the addition makes sense for the validity of the matrix.
   * For example, a right-to-left addition could corrupt the computation of the barcode if done blindly.
   * So should be used with care.
   * 
   * @param coefficient Value to multiply.
   * @param sourcePosition @ref PosIdx index of the source column.
   * @param targetPosition @ref PosIdx index of the target column.
   */
  void multiply_source_and_add_to(const Field_element& coefficient, 
                                  Pos_index sourcePosition,
                                  Pos_index targetPosition);

  /**
   * @brief Indicates if the entry at given coordinates has value zero.
   * 
   * @param position @ref PosIdx index of the cell corresponding to the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @return true If the entry has value zero.
   * @return false Otherwise.
   */
  bool is_zero_entry(Pos_index position, ID_index rowIndex) const;
  /**
   * @brief Indicates if the column at given index has value zero.
   *
   * Note that this method should always return false, as a valid @ref chainmatrix "chain matrix" never has
   * empty columns.
   * 
   * @param position @ref PosIdx index of the cell corresponding to the column.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(Pos_index position);

  /**
   * @brief Returns the @ref PosIdx index of the column which has the given @ref rowindex "row index" as pivot.
   * Assumes that the pivot exists.
   * 
   * @param cellIndex @ref rowindex "Row index" of the pivot.
   * @return @ref PosIdx index of the column with the given pivot.
   */
  Pos_index get_column_with_pivot(ID_index cellIndex) const;  // assumes that pivot exists
  /**
   * @brief Returns the @ref rowindex "row index" of the pivot of the given column.
   * 
   * @param position @ref PosIdx index of the cell corresponding to the column.
   * @return The @ref rowindex "row index" of the pivot.
   */
  ID_index get_pivot(Pos_index position);

  /**
   * @brief Resets the matrix to an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  void reset(Column_settings* colSettings) {
    matrix_.reset(colSettings);
    positionToIndex_.clear();
    nextPosition_ = 0;
    nextIndex_ = 0;
  }

  // void set_operators(Field_operators* operators) { matrix_.set_operators(operators); }

  /**
   * @brief Assign operator.
   */
  Position_to_index_overlay& operator=(const Position_to_index_overlay& other);
  /**
   * @brief Swap operator.
   */
  friend void swap(Position_to_index_overlay& matrix1, Position_to_index_overlay& matrix2) {
    swap(matrix1.matrix_, matrix2.matrix_);
    matrix1.positionToIndex_.swap(matrix2.positionToIndex_);
    std::swap(matrix1.nextPosition_, matrix2.nextPosition_);
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
   * @return A reference to the barcode. The barcode is a vector of @ref Matrix::Bar. A bar stores three informations:
   * the @ref PosIdx birth index, the @ref PosIdx death index and the dimension of the bar.
   */
  const Barcode& get_current_barcode() const;
  
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
  
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true.
   * Does the same than @ref vine_swap, but assumes that the swap is non trivial and
   * therefore skips a part of the case study.
   * 
   * @param position @ref PosIdx index of the first cell to swap. The second one has to be at `position + 1`.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap_with_z_eq_1_case(Pos_index position);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true.
   * Does a vine swap between two cells which are consecutive in the filtration. Roughly, if \f$ F \f$ is the current
   * filtration represented by the matrix, the method modifies the matrix such that the new state corresponds to 
   * a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but with the two cells at position `position`
   * and `position + 1` swapped. Of course, the two cells should not have a face/coface relation which each other ;
   * \f$ F' \f$ has to be a valid filtration.
   * See @cite vineyards for more information about vine and vineyards.
   * 
   * @param position @ref PosIdx index of the first cell to swap. The second one has to be at `position + 1`.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap(Pos_index position);

 private:
  Underlying_matrix matrix_;            /**< Interfaced matrix. */
  std::vector<Index> positionToIndex_;  /**< Map from @ref PosIdx index to @ref MatIdx index. */
  Pos_index nextPosition_;              /**< Next unused position. */
  Index nextIndex_;                     /**< Next unused index. */
};

template <class Underlying_matrix, class Master_matrix>
inline Position_to_index_overlay<Underlying_matrix, Master_matrix>::Position_to_index_overlay(
    Column_settings* colSettings)
    : matrix_(colSettings), nextPosition_(0), nextIndex_(0) 
{}

template <class Underlying_matrix, class Master_matrix>
template <class Boundary_range>
inline Position_to_index_overlay<Underlying_matrix, Master_matrix>::Position_to_index_overlay(
    const std::vector<Boundary_range>& orderedBoundaries, Column_settings* colSettings)
    : matrix_(orderedBoundaries, colSettings),
      positionToIndex_(orderedBoundaries.size()),
      nextPosition_(orderedBoundaries.size()),
      nextIndex_(orderedBoundaries.size()) 
{
  for (Index i = 0; i < orderedBoundaries.size(); i++) {
    positionToIndex_[i] = i;
  }
}

template <class Underlying_matrix, class Master_matrix>
inline Position_to_index_overlay<Underlying_matrix, Master_matrix>::Position_to_index_overlay(
    unsigned int numberOfColumns, Column_settings* colSettings)
    : matrix_(numberOfColumns, colSettings),
      positionToIndex_(numberOfColumns),
      nextPosition_(0),
      nextIndex_(0) 
{}

template <class Underlying_matrix, class Master_matrix>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Position_to_index_overlay<Underlying_matrix, Master_matrix>::Position_to_index_overlay(
    Column_settings* colSettings, 
    const BirthComparatorFunction& birthComparator, 
    const DeathComparatorFunction& deathComparator)
    : matrix_(colSettings, birthComparator, deathComparator), nextPosition_(0), nextIndex_(0) 
{}

template <class Underlying_matrix, class Master_matrix>
template <typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_range>
inline Position_to_index_overlay<Underlying_matrix, Master_matrix>::Position_to_index_overlay(
    const std::vector<Boundary_range>& orderedBoundaries, 
    Column_settings* colSettings,
    const BirthComparatorFunction& birthComparator, 
    const DeathComparatorFunction& deathComparator)
    : matrix_(orderedBoundaries, colSettings, birthComparator, deathComparator),
      positionToIndex_(orderedBoundaries.size()),
      nextPosition_(orderedBoundaries.size()),
      nextIndex_(orderedBoundaries.size()) 
{
  for (Index i = 0; i < orderedBoundaries.size(); i++) {
    positionToIndex_[i] = i;
  }
}

template <class Underlying_matrix, class Master_matrix>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Position_to_index_overlay<Underlying_matrix, Master_matrix>::Position_to_index_overlay(
    unsigned int numberOfColumns, 
    Column_settings* colSettings,
    const BirthComparatorFunction& birthComparator, 
    const DeathComparatorFunction& deathComparator)
    : matrix_(numberOfColumns, colSettings, birthComparator, deathComparator),
      positionToIndex_(numberOfColumns),
      nextPosition_(0),
      nextIndex_(0) 
{}

template <class Underlying_matrix, class Master_matrix>
inline Position_to_index_overlay<Underlying_matrix, Master_matrix>::Position_to_index_overlay(
    const Position_to_index_overlay& matrixToCopy, Column_settings* colSettings)
    : matrix_(matrixToCopy.matrix_, colSettings),
      positionToIndex_(matrixToCopy.positionToIndex_),
      nextPosition_(matrixToCopy.nextPosition_),
      nextIndex_(matrixToCopy.nextIndex_) 
{}

template <class Underlying_matrix, class Master_matrix>
inline Position_to_index_overlay<Underlying_matrix, Master_matrix>::Position_to_index_overlay(
    Position_to_index_overlay&& other) noexcept
    : matrix_(std::move(other.matrix_)),
      positionToIndex_(std::move(other.positionToIndex_)),
      nextPosition_(std::exchange(other.nextPosition_, 0)),
      nextIndex_(std::exchange(other.nextIndex_, 0)) 
{}

template <class Underlying_matrix, class Master_matrix>
template <class Boundary_range>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::insert_boundary(const Boundary_range& boundary,
                                                                                        Dimension dim) 
{
  if (positionToIndex_.size() <= nextPosition_) {
    positionToIndex_.resize(nextPosition_ * 2 + 1);
  }

  positionToIndex_[nextPosition_++] = nextIndex_++;

  matrix_.insert_boundary(boundary, dim);
}

template <class Underlying_matrix, class Master_matrix>
template <class Boundary_range>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::insert_boundary(ID_index cellIndex,
                                                                                        const Boundary_range& boundary,
                                                                                        Dimension dim) 
{
  if (positionToIndex_.size() <= nextPosition_) {
    positionToIndex_.resize(nextPosition_ * 2 + 1);
  }

  positionToIndex_[nextPosition_++] = nextIndex_++;

  matrix_.insert_boundary(cellIndex, boundary, dim);
}

template <class Underlying_matrix, class Master_matrix>
inline typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Column&
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_column(Pos_index position) 
{
  return matrix_.get_column(positionToIndex_[position]);
}

template <class Underlying_matrix, class Master_matrix>
inline const typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Column&
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_column(Pos_index position) const 
{
  return matrix_.get_column(positionToIndex_[position]);
}

template <class Underlying_matrix, class Master_matrix>
inline typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Row&
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_row(ID_index rowIndex) 
{
  return matrix_.get_row(rowIndex);
}

template <class Underlying_matrix, class Master_matrix>
inline const typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Row&
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_row(ID_index rowIndex) const 
{
  return matrix_.get_row(rowIndex);
}

template <class Underlying_matrix, class Master_matrix>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::erase_empty_row(ID_index rowIndex) 
{
  return matrix_.erase_empty_row(rowIndex);
}

template <class Underlying_matrix, class Master_matrix>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::remove_maximal_cell(Pos_index position) 
{
  --nextPosition_;

  ID_index pivot = matrix_.get_pivot(positionToIndex_[position]);
  std::vector<Index> columnsToSwap(nextPosition_ - position);

  if (nextPosition_ != position) {
    positionToIndex_[position] = positionToIndex_[position + 1];
    for (Pos_index p = position + 1; p < nextPosition_; ++p) {
      columnsToSwap[p - position - 1] = positionToIndex_[p];
      positionToIndex_[p] = positionToIndex_[p + 1];
    }
    columnsToSwap.back() = positionToIndex_[nextPosition_];
  }

  matrix_.remove_maximal_cell(pivot, columnsToSwap);
}

template <class Underlying_matrix, class Master_matrix>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::remove_last() 
{
  --nextPosition_;
  if constexpr (Master_matrix::Option_list::has_vine_update) {
    std::vector<Index> columnsToSwap;
    matrix_.remove_maximal_cell(matrix_.get_pivot(positionToIndex_[nextPosition_]), columnsToSwap);
  } else {
    matrix_.remove_last();  // linear with vine updates, so it is better to use remove_maximal_cell
  }
}

template <class Underlying_matrix, class Master_matrix>
inline typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Dimension
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_max_dimension() const 
{
  return matrix_.get_max_dimension();
}

template <class Underlying_matrix, class Master_matrix>
inline typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Index
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_number_of_columns() const 
{
  return matrix_.get_number_of_columns();
}

template <class Underlying_matrix, class Master_matrix>
inline typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Dimension
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_column_dimension(Pos_index position) const 
{
  return matrix_.get_column_dimension(positionToIndex_[position]);
}

template <class Underlying_matrix, class Master_matrix>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::add_to(Pos_index sourcePosition,
                                                                               Pos_index targetPosition) 
{
  return matrix_.add_to(positionToIndex_[sourcePosition], positionToIndex_[targetPosition]);
}

template <class Underlying_matrix, class Master_matrix>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::multiply_target_and_add_to(
    Pos_index sourcePosition, const Field_element& coefficient, Pos_index targetPosition) 
{
  return matrix_.multiply_target_and_add_to(positionToIndex_[sourcePosition], 
                                            coefficient,
                                            positionToIndex_[targetPosition]);
}

template <class Underlying_matrix, class Master_matrix>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::multiply_source_and_add_to(
    const Field_element& coefficient, Pos_index sourcePosition, Pos_index targetPosition) 
{
  return matrix_.multiply_source_and_add_to(coefficient, 
                                            positionToIndex_[sourcePosition],
                                            positionToIndex_[targetPosition]);
}

template <class Underlying_matrix, class Master_matrix>
inline bool Position_to_index_overlay<Underlying_matrix, Master_matrix>::is_zero_entry(Pos_index position,
                                                                                     ID_index rowIndex) const 
{
  return matrix_.is_zero_entry(positionToIndex_[position], rowIndex);
}

template <class Underlying_matrix, class Master_matrix>
inline bool Position_to_index_overlay<Underlying_matrix, Master_matrix>::is_zero_column(Pos_index position) 
{
  return matrix_.is_zero_column(positionToIndex_[position]);
}

template <class Underlying_matrix, class Master_matrix>
inline typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Pos_index
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_column_with_pivot(ID_index cellIndex) const 
{
  Index id = matrix_.get_column_with_pivot(cellIndex);
  Pos_index i = 0;
  while (positionToIndex_[i] != id) ++i;
  return i;
}

template <class Underlying_matrix, class Master_matrix>
inline typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::ID_index
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_pivot(Pos_index position) 
{
  return matrix_.get_pivot(positionToIndex_[position]);
}

template <class Underlying_matrix, class Master_matrix>
inline Position_to_index_overlay<Underlying_matrix, Master_matrix>&
Position_to_index_overlay<Underlying_matrix, Master_matrix>::operator=(const Position_to_index_overlay& other) 
{
  matrix_ = other.matrix_;
  positionToIndex_ = other.positionToIndex_;
  nextPosition_ = other.nextPosition_;
  nextIndex_ = other.nextIndex_;

  return *this;
}

template <class Underlying_matrix, class Master_matrix>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::print() 
{
  return matrix_.print();
}

template <class Underlying_matrix, class Master_matrix>
inline const typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Barcode&
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_current_barcode() const 
{
  return matrix_.get_current_barcode();
}

template <class Underlying_matrix, class Master_matrix>
inline void Position_to_index_overlay<Underlying_matrix, Master_matrix>::update_representative_cycles() 
{
  matrix_.update_representative_cycles();
}

template <class Underlying_matrix, class Master_matrix>
inline const std::vector<typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Cycle>&
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_representative_cycles() 
{
  return matrix_.get_representative_cycles();
}

template <class Underlying_matrix, class Master_matrix>
inline const typename Position_to_index_overlay<Underlying_matrix, Master_matrix>::Cycle&
Position_to_index_overlay<Underlying_matrix, Master_matrix>::get_representative_cycle(const Bar& bar) 
{
  return matrix_.get_representative_cycle(bar);
}

template <class Underlying_matrix, class Master_matrix>
inline bool Position_to_index_overlay<Underlying_matrix, Master_matrix>::vine_swap_with_z_eq_1_case(Pos_index position) 
{
  Index next = matrix_.vine_swap_with_z_eq_1_case(positionToIndex_[position], positionToIndex_[position + 1]);
  if (next == positionToIndex_[position]) {
    std::swap(positionToIndex_[position], positionToIndex_[position + 1]);
    return true;
  }

  return false;
}

template <class Underlying_matrix, class Master_matrix>
inline bool Position_to_index_overlay<Underlying_matrix, Master_matrix>::vine_swap(Pos_index position) 
{
  Index next = matrix_.vine_swap(positionToIndex_[position], positionToIndex_[position + 1]);
  if (next == positionToIndex_[position]) {
    std::swap(positionToIndex_[position], positionToIndex_[position + 1]);
    return true;
  }

  return false;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_POS_TO_ID_TRANSLATION_H
