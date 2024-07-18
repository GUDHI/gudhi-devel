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
 * @file overlay_posidx_to_matidx.h
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
 * @class Position_to_index_overlay overlay_posidx_to_matidx.h gudhi/Persistence_matrix/overlay_posidx_to_matidx.h
 * @ingroup persistence_matrix
 *
 * @brief Overlay for @ref chainmatrix "chain matrices" replacing all input and output @ref MatIdx indices of the
 * original methods with @ref PosIdx indices. The overlay is useless for @ref boundarymatrix "boundary matrices"
 * as @ref MatIdx == @ref PosIdx for them.
 * 
 * @tparam %Matrix_type Matrix type taking the overlay.
 * @tparam Master_matrix_type An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Matrix_type, class Master_matrix_type>
class Position_to_index_overlay 
{
 public:
  using index = typename Master_matrix_type::index;                       /**< @ref MatIdx index type. */
  using id_index = typename Master_matrix_type::id_index;                 /**< @ref IDIdx index type. */
  using pos_index = typename Master_matrix_type::pos_index;               /**< @ref PosIdx index type. */
  using dimension_type = typename Master_matrix_type::dimension_type;     /**< Dimension value type. */
  /**
   * @brief Field operators class. Necessary only if @ref PersistenceMatrixOptions::is_z2 is false.
   */
  using Field_operators = typename Master_matrix_type::Field_operators;
  using Field_element_type = typename Master_matrix_type::element_type;   /**< Type of an field element. */
  using boundary_type = typename Master_matrix_type::boundary_type;       /**< Type of an input column. */
  using Column_type = typename Master_matrix_type::Column_type;           /**< Column type. */
  using Row_type = typename Master_matrix_type::Row_type;                 /**< Row type,
                                                                               only necessary with row access option. */
  using bar_type = typename Master_matrix_type::Bar;                      /**< Bar type. */
  using barcode_type = typename Master_matrix_type::barcode_type;         /**< Barcode type. */
  using cycle_type = typename Master_matrix_type::cycle_type;             /**< Cycle type. */
  using cell_rep_type = typename Master_matrix_type::cell_rep_type;       /**< %Cell content representative. */
  using Cell_constructor = typename Master_matrix_type::Cell_constructor; /**< Factory of @ref Cell classes. */
  using Column_settings = typename Master_matrix_type::Column_settings;   /**< Structure giving access to the columns to
                                                                               necessary external classes. */

  /**
   * @brief Constructs an empty matrix.
   * 
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  Position_to_index_overlay(Column_settings* colSettings);
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
   * matrix with the @ref Position_to_index_overlay(unsigned int, Column_settings*)
   * constructor preferably).
   * @param colSettings Pointer to an existing setting structure for the columns. The structure should contain all
   * the necessary external classes specifically necessary for the choosen column type, such as custom allocators.
   */
  template <class Boundary_type = boundary_type>
  Position_to_index_overlay(const std::vector<Boundary_type>& orderedBoundaries, 
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
   * @tparam BirthComparatorFunction Type of the birth comparator: (@ref pos_index, @ref pos_index) -> bool
   * @tparam DeathComparatorFunction Type of the death comparator: (@ref pos_index, @ref pos_index) -> bool
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
   * Constructs a new matrix from the given ranges of @ref Matrix::cell_rep_type. Each range corresponds to a column 
   * (the order of the ranges are preserved). The content of the ranges is assumed to be sorted by increasing IDs.
   * The IDs of the simplices are also assumed to be consecutive, ordered by filtration value, starting with 0. 
   *
   * @warning If @ref PersistenceMatrixOptions::has_vine_update is false, the comparators are not used.
   * And if @ref PersistenceMatrixOptions::has_vine_update is true, but
   * @ref PersistenceMatrixOptions::has_column_pairings is also true, the comparators are ignored and
   * the current barcode is used to compare birth and deaths. Therefore it is useless to provide them in those cases.
   * 
   * @tparam BirthComparatorFunction Type of the birth comparator: (@ref pos_index, @ref pos_index) -> bool
   * @tparam DeathComparatorFunction Type of the death comparator: (@ref pos_index, @ref pos_index) -> bool
   * @tparam Boundary_type  Range type for @ref Matrix::cell_rep_type ranges.
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
  template <typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_type>
  Position_to_index_overlay(const std::vector<Boundary_type>& orderedBoundaries, 
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
   * @tparam BirthComparatorFunction Type of the birth comparator: (@ref pos_index, @ref pos_index) -> bool
   * @tparam DeathComparatorFunction Type of the death comparator: (@ref pos_index, @ref pos_index) -> bool
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
   * It also assumes that the faces in the given boundary are identified by their relative position in the filtration, 
   * starting at 0. If it is not the case, use the other
   * @ref insert_boundary(id_index, const Boundary_type&, dimension_type) "insert_boundary" instead by indicating the
   * face ID used in the boundaries when the face is inserted.
   *
   * Different to the constructor, the boundaries do not have to come from a simplicial complex, but also from
   * a more general cell complex. This includes cubical complexes or Morse complexes for example.
   *
   * When inserted, the given boundary is reduced and from the reduction process, the column is deduced in the form of:
   * `IDIdx + linear combination of older column IDIdxs`. If the barcode is stored, it will be updated.
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
   * @p faceID values of precedent calls of the method for the corresponding faces and should be ordered in 
   * increasing order.
   * @param dim Dimension of the face whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   */
  template <class Boundary_type = boundary_type>
  void insert_boundary(id_index faceIndex, const Boundary_type& boundary, dimension_type dim = -1);
  /**
   * @brief Returns the column at the given @ref PosIdx index.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param position @ref PosIdx index of the column to return.
   * @return Reference to the column.
   */
  Column_type& get_column(pos_index position);
  /**
   * @brief Returns the column at the given @ref PosIdx index.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param position @ref PosIdx index of the column to return.
   * @return Const reference to the column.
   */
  const Column_type& get_column(pos_index position) const;
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access is true.
   * Returns the row at the given @ref rowindex "row index".
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return.
   * @return Reference to the row.
   */
  Row_type& get_row(id_index rowIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access is true.
   * Returns the row at the given @ref rowindex "row index".
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * 
   * @param rowIndex @ref rowindex "Row index" of the row to return.
   * @return Const reference to the row.
   */
  const Row_type& get_row(id_index rowIndex) const;
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access and
   * @ref PersistenceMatrixOptions::has_removable_rows are true.
   * Assumes that the row is empty and removes it. 
   *
   * @warning The removed rows are always assumed to be empty. If it is not the case, the deleted row cells are not
   * removed from their columns. And in the case of intrusive rows, this will generate a segmentation fault when 
   * the column cells are destroyed later. The row access is just meant as a "read only" access to the rows and the
   * @ref erase_empty_row method just as a way to specify that a row is empty and can therefore be removed from
   * dictionaries. This allows to avoid testing the emptiness of a row at each column cell removal, what can be
   * quite frequent. 
   * 
   * @param rowIndex @ref rowindex "Row index" of the empty row to remove.
   */
  void erase_empty_row(id_index rowIndex);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns,
   * @ref PersistenceMatrixOptions::has_vine_update and @ref PersistenceMatrixOptions::has_map_column_container
   * are true.
   * Assumes that the face is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., the matrix is still a compatible bases of the chain complex in the sense of @cite zigzag).
   * The maximality of the face is not verified.
   * Also updates the barcode if it was computed.
   *
   * See also @ref remove_last.
   * 
   * @param position @ref PosIdx index of the face to remove.
   */
  void remove_maximal_face(pos_index position);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_removable_columns is true and,
   * if @ref PersistenceMatrixOptions::has_map_column_container is true or
   * @ref PersistenceMatrixOptions::has_vine_update is false.
   * Removes the last face in the filtration from the matrix and updates the barcode if it is stored.
   *
   * See also @ref remove_maximal_face.
   */
  void remove_last();

  /**
   * @brief Returns the maximal dimension of a face stored in the matrix. Only available 
   * if @ref PersistenceMatrixOptions::has_matrix_maximal_dimension_access is true.
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
   * @brief Returns the dimension of the given face.
   * 
   * @param position @ref PosIdx index of the face.
   * @return Dimension of the face.
   */
  dimension_type get_column_dimension(pos_index position) const;

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
  void add_to(pos_index sourcePosition, pos_index targetPosition);
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
  void multiply_target_and_add_to(pos_index sourcePosition, 
                                  const Field_element_type& coefficient,
                                  pos_index targetPosition);
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
  void multiply_source_and_add_to(const Field_element_type& coefficient, 
                                  pos_index sourcePosition,
                                  pos_index targetPosition);

  /**
   * @brief Indicates if the cell at given coordinates has value zero.
   * 
   * @param position @ref PosIdx index of the face corresponding to the column of the cell.
   * @param rowIndex @ref rowindex "Row index" of the row of the cell.
   * @return true If the cell has value zero.
   * @return false Otherwise.
   */
  bool is_zero_cell(pos_index position, id_index rowIndex) const;
  /**
   * @brief Indicates if the column at given index has value zero.
   *
   * Note that this method should always return false, as a valid @ref chainmatrix "chain matrix" never has
   * empty columns.
   * 
   * @param position @ref PosIdx index of the face corresponding to the column.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(pos_index position);

  /**
   * @brief Returns the @ref PosIdx index of the column which has the given @ref rowindex "row index" as pivot.
   * Assumes that the pivot exists.
   * 
   * @param faceIndex @ref rowindex "Row index" of the pivot.
   * @return @ref PosIdx index of the column with the given pivot.
   */
  pos_index get_column_with_pivot(id_index faceIndex) const;  // assumes that pivot exists
  /**
   * @brief Returns the @ref rowindex "row index" of the pivot of the given column.
   * 
   * @param position @ref PosIdx index of the face corresponding to the column.
   * @return The @ref rowindex "row index" of the pivot.
   */
  id_index get_pivot(pos_index position);

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
  const barcode_type& get_current_barcode() const;
  
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
  const std::vector<cycle_type>& get_representative_cycles();
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true.
   * Returns the cycle representing the given bar.
   * 
   * @param bar A bar from the current barcode.
   * @return A const reference to the cycle representing @p bar.
   */
  const cycle_type& get_representative_cycle(const bar_type& bar);
  
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true.
   * Does the same than @ref vine_swap, but assumes that the swap is non trivial and
   * therefore skips a part of the case study.
   * 
   * @param position @ref PosIdx index of the first face to swap. The second one has to be at `position + 1`.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap_with_z_eq_1_case(pos_index position);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true.
   * Does a vine swap between two faces which are consecutive in the filtration. Roughly, if \f$ F \f$ is the current
   * filtration represented by the matrix, the method modifies the matrix such that the new state corresponds to 
   * a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but with the two faces at position `position`
   * and `position + 1` swapped. Of course, the two faces should not have a face/coface relation which each other ;
   * \f$ F' \f$ has to be a valid filtration.
   * See @cite vineyards for more information about vine and vineyards.
   * 
   * @param position @ref PosIdx index of the first face to swap. The second one has to be at `position + 1`.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap(pos_index position);

 private:
  Matrix_type matrix_;                  /**< Interfaced matrix. */
  std::vector<index> positionToIndex_;  /**< Map from @ref PosIdx index to @ref MatIdx index. */
  pos_index nextPosition_;              /**< Next unused position. */
  index nextIndex_;                     /**< Next unused index. */
};

template <class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type, Master_matrix_type>::Position_to_index_overlay(
    Column_settings* colSettings)
    : matrix_(colSettings), nextPosition_(0), nextIndex_(0) 
{}

template <class Matrix_type, class Master_matrix_type>
template <class Boundary_type>
inline Position_to_index_overlay<Matrix_type, Master_matrix_type>::Position_to_index_overlay(
    const std::vector<Boundary_type>& orderedBoundaries, Column_settings* colSettings)
    : matrix_(orderedBoundaries, colSettings),
      positionToIndex_(orderedBoundaries.size()),
      nextPosition_(orderedBoundaries.size()),
      nextIndex_(orderedBoundaries.size()) 
{
  for (index i = 0; i < orderedBoundaries.size(); i++) {
    positionToIndex_[i] = i;
  }
}

template <class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type, Master_matrix_type>::Position_to_index_overlay(
    unsigned int numberOfColumns, Column_settings* colSettings)
    : matrix_(numberOfColumns, colSettings),
      positionToIndex_(numberOfColumns),
      nextPosition_(0),
      nextIndex_(0) 
{}

template <class Matrix_type, class Master_matrix_type>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Position_to_index_overlay<Matrix_type, Master_matrix_type>::Position_to_index_overlay(
    Column_settings* colSettings, 
    const BirthComparatorFunction& birthComparator, 
    const DeathComparatorFunction& deathComparator)
    : matrix_(colSettings, birthComparator, deathComparator), nextPosition_(0), nextIndex_(0) 
{}

template <class Matrix_type, class Master_matrix_type>
template <typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_type>
inline Position_to_index_overlay<Matrix_type, Master_matrix_type>::Position_to_index_overlay(
    const std::vector<Boundary_type>& orderedBoundaries, 
    Column_settings* colSettings,
    const BirthComparatorFunction& birthComparator, 
    const DeathComparatorFunction& deathComparator)
    : matrix_(orderedBoundaries, colSettings, birthComparator, deathComparator),
      positionToIndex_(orderedBoundaries.size()),
      nextPosition_(orderedBoundaries.size()),
      nextIndex_(orderedBoundaries.size()) 
{
  for (index i = 0; i < orderedBoundaries.size(); i++) {
    positionToIndex_[i] = i;
  }
}

template <class Matrix_type, class Master_matrix_type>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Position_to_index_overlay<Matrix_type, Master_matrix_type>::Position_to_index_overlay(
    unsigned int numberOfColumns, 
    Column_settings* colSettings,
    const BirthComparatorFunction& birthComparator, 
    const DeathComparatorFunction& deathComparator)
    : matrix_(numberOfColumns, colSettings, birthComparator, deathComparator),
      positionToIndex_(numberOfColumns),
      nextPosition_(0),
      nextIndex_(0) 
{}

template <class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type, Master_matrix_type>::Position_to_index_overlay(
    const Position_to_index_overlay& matrixToCopy, Column_settings* colSettings)
    : matrix_(matrixToCopy.matrix_, colSettings),
      positionToIndex_(matrixToCopy.positionToIndex_),
      nextPosition_(matrixToCopy.nextPosition_),
      nextIndex_(matrixToCopy.nextIndex_) 
{}

template <class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type, Master_matrix_type>::Position_to_index_overlay(
    Position_to_index_overlay&& other) noexcept
    : matrix_(std::move(other.matrix_)),
      positionToIndex_(std::move(other.positionToIndex_)),
      nextPosition_(std::exchange(other.nextPosition_, 0)),
      nextIndex_(std::exchange(other.nextIndex_, 0)) 
{}

template <class Matrix_type, class Master_matrix_type>
template <class Boundary_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::insert_boundary(const Boundary_type& boundary,
                                                                                        dimension_type dim) 
{
  if (positionToIndex_.size() <= nextPosition_) {
    positionToIndex_.resize(nextPosition_ * 2 + 1);
  }

  positionToIndex_[nextPosition_++] = nextIndex_++;

  matrix_.insert_boundary(boundary, dim);
}

template <class Matrix_type, class Master_matrix_type>
template <class Boundary_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::insert_boundary(id_index faceIndex,
                                                                                        const Boundary_type& boundary,
                                                                                        dimension_type dim) 
{
  if (positionToIndex_.size() <= nextPosition_) {
    positionToIndex_.resize(nextPosition_ * 2 + 1);
  }

  positionToIndex_[nextPosition_++] = nextIndex_++;

  matrix_.insert_boundary(faceIndex, boundary, dim);
}

template <class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::Column_type&
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_column(pos_index position) 
{
  return matrix_.get_column(positionToIndex_[position]);
}

template <class Matrix_type, class Master_matrix_type>
inline const typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::Column_type&
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_column(pos_index position) const 
{
  return matrix_.get_column(positionToIndex_[position]);
}

template <class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::Row_type&
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_row(id_index rowIndex) 
{
  return matrix_.get_row(rowIndex);
}

template <class Matrix_type, class Master_matrix_type>
inline const typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::Row_type&
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_row(id_index rowIndex) const 
{
  return matrix_.get_row(rowIndex);
}

template <class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::erase_empty_row(id_index rowIndex) 
{
  return matrix_.erase_empty_row(rowIndex);
}

template <class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::remove_maximal_face(pos_index position) 
{
  --nextPosition_;

  id_index pivot = matrix_.get_pivot(positionToIndex_[position]);
  std::vector<index> columnsToSwap(nextPosition_ - position);

  if (nextPosition_ != position) {
    positionToIndex_[position] = positionToIndex_[position + 1];
    for (pos_index p = position + 1; p < nextPosition_; ++p) {
      columnsToSwap[p - position - 1] = positionToIndex_[p];
      positionToIndex_[p] = positionToIndex_[p + 1];
    }
    columnsToSwap.back() = positionToIndex_[nextPosition_];
  }

  matrix_.remove_maximal_face(pivot, columnsToSwap);
}

template <class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::remove_last() 
{
  --nextPosition_;
  if constexpr (Master_matrix_type::Option_list::has_vine_update) {
    std::vector<index> columnsToSwap;
    matrix_.remove_maximal_face(matrix_.get_pivot(positionToIndex_[nextPosition_]), columnsToSwap);
  } else {
    matrix_.remove_last();  // linear with vine updates, so it is better to use remove_maximal_face
  }
}

template <class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::dimension_type
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_max_dimension() const 
{
  return matrix_.get_max_dimension();
}

template <class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::index
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_number_of_columns() const 
{
  return matrix_.get_number_of_columns();
}

template <class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::dimension_type
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_column_dimension(pos_index position) const 
{
  return matrix_.get_column_dimension(positionToIndex_[position]);
}

template <class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::add_to(pos_index sourcePosition,
                                                                               pos_index targetPosition) 
{
  return matrix_.add_to(positionToIndex_[sourcePosition], positionToIndex_[targetPosition]);
}

template <class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::multiply_target_and_add_to(
    pos_index sourcePosition, const Field_element_type& coefficient, pos_index targetPosition) 
{
  return matrix_.multiply_target_and_add_to(positionToIndex_[sourcePosition], 
                                            coefficient,
                                            positionToIndex_[targetPosition]);
}

template <class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::multiply_source_and_add_to(
    const Field_element_type& coefficient, pos_index sourcePosition, pos_index targetPosition) 
{
  return matrix_.multiply_source_and_add_to(coefficient, 
                                            positionToIndex_[sourcePosition],
                                            positionToIndex_[targetPosition]);
}

template <class Matrix_type, class Master_matrix_type>
inline bool Position_to_index_overlay<Matrix_type, Master_matrix_type>::is_zero_cell(pos_index position,
                                                                                     id_index rowIndex) const 
{
  return matrix_.is_zero_cell(positionToIndex_[position], rowIndex);
}

template <class Matrix_type, class Master_matrix_type>
inline bool Position_to_index_overlay<Matrix_type, Master_matrix_type>::is_zero_column(pos_index position) 
{
  return matrix_.is_zero_column(positionToIndex_[position]);
}

template <class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::pos_index
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_column_with_pivot(id_index faceIndex) const 
{
  index id = matrix_.get_column_with_pivot(faceIndex);
  pos_index i = 0;
  while (positionToIndex_[i] != id) ++i;
  return i;
}

template <class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::id_index
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_pivot(pos_index position) 
{
  return matrix_.get_pivot(positionToIndex_[position]);
}

template <class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type, Master_matrix_type>&
Position_to_index_overlay<Matrix_type, Master_matrix_type>::operator=(const Position_to_index_overlay& other) 
{
  matrix_ = other.matrix_;
  positionToIndex_ = other.positionToIndex_;
  nextPosition_ = other.nextPosition_;
  nextIndex_ = other.nextIndex_;

  return *this;
}

template <class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::print() 
{
  return matrix_.print();
}

template <class Matrix_type, class Master_matrix_type>
inline const typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::barcode_type&
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_current_barcode() const 
{
  return matrix_.get_current_barcode();
}

template <class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type, Master_matrix_type>::update_representative_cycles() 
{
  matrix_.update_representative_cycles();
}

template <class Matrix_type, class Master_matrix_type>
inline const std::vector<typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::cycle_type>&
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_representative_cycles() 
{
  return matrix_.get_representative_cycles();
}

template <class Matrix_type, class Master_matrix_type>
inline const typename Position_to_index_overlay<Matrix_type, Master_matrix_type>::cycle_type&
Position_to_index_overlay<Matrix_type, Master_matrix_type>::get_representative_cycle(const bar_type& bar) 
{
  return matrix_.get_representative_cycle(bar);
}

template <class Matrix_type, class Master_matrix_type>
inline bool Position_to_index_overlay<Matrix_type, Master_matrix_type>::vine_swap_with_z_eq_1_case(pos_index position) 
{
  index next = matrix_.vine_swap_with_z_eq_1_case(positionToIndex_[position], positionToIndex_[position + 1]);
  if (next == positionToIndex_[position]) {
    std::swap(positionToIndex_[position], positionToIndex_[position + 1]);
    return true;
  }

  return false;
}

template <class Matrix_type, class Master_matrix_type>
inline bool Position_to_index_overlay<Matrix_type, Master_matrix_type>::vine_swap(pos_index position) 
{
  index next = matrix_.vine_swap(positionToIndex_[position], positionToIndex_[position + 1]);
  if (next == positionToIndex_[position]) {
    std::swap(positionToIndex_[position], positionToIndex_[position + 1]);
    return true;
  }

  return false;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_POS_TO_ID_TRANSLATION_H
