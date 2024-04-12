/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file PersistenceMatrixOptions.h
 * @author Hannah Schreiber
 * @brief Contains the concept for the matrix options.
 */

/// Gudhi namespace.
namespace Gudhi {
/// Persistence matrix namespace.
namespace persistence_matrix {

/** 
 * @ingroup persistence_matrix
 *
 * @brief Concept of the template parameter for the class @ref Matrix.
 *
 * An implementation of this concept is @ref Default_options.
 * If you want to provide your own, it is recommended that you derive from it and override some parts instead of
 * writing a class from scratch.
 */
struct PersistenceMatrixOptions 
{
  /**
   * @brief Field operators. Has to follow the @ref FieldOperators concept.
   * The type will not be used if @ref is_z2 is set to true, so it can be set to anything.
   */
  using Field_coeff_operators = unspecified;
  /**
   * @brief Type for the dimension. Has to be an integer type.
   * If unsigned, the maximal value of the type should not be attained during a run.
   */
  using dimension_type = unspecified;
  /**
   * @brief Type for the different indexation types and should be able to contain the maximal number of columns
   * of the matrix during a run. Has to be an integer type.
   * If unsigned, the maximal value of the type should not be attained during a run.
   */
  using index_type = unspecified;

  /**
   * @brief If true, indicates that the values contained in the matrix are in \f$ Z_2 \f$ and can therefore
   * be treated like booleans. If set to false, the values are assumed to be in the field \f$ Z_p \f$ for
   * some prime \f$ p \f$ given by @ref Field_coeff_operators. It is highly recommended to set the variable to true,
   * if \f$ p = 2 \f$.
   */
  static const bool is_z2;
  /**
   * @brief Specifies the desired column type. All possible column types are described in @ref Column_types.
   */
  static const Column_types column_type;
  /**
   * @brief Specifies the desired indexation scheme to access the methods of the matrix.
   * See @ref mp_indexation "matrix description" and @ref Column_indexation_types for more details about the meaning
   * of the indexation types.
   */
  static const Column_indexation_types column_indexation_type;

  /**
   * @brief Only enabled for @ref basematrix "base matrices" (i.e., none of the following is true:
   * @ref has_column_pairings, @ref has_vine_update, @ref can_retrieve_representative_cycles), is ignored otherwise.
   * If set to true, two identical columns in the matrix are not explicitely stored separately but are represented 
   * by a same column.
   *
   * Note that some methods of the @ref basematrix "base matrix" are not available when true:
   * - @ref Matrix::insert_column(const Container_type&, index_type) "insert_column(const Container_type&, index)",
   * - @ref Matrix::zero_column(index_type) "zero_column(index)",
   * - @ref Matrix::zero_cell(index_type, index_type) "zero_cell(index, id_index)",
   * - @ref Matrix::swap_columns(index_type, index_type) "swap_columns(index, index)",
   * - @ref Matrix::swap_rows(index_type, index_type) "swap_rows(index, index)",
   * - @ref Matrix::remove_column(index_type) "remove_column(index)",
   * - @ref Matrix::remove_last "remove_last()".
   */
  static const bool has_column_compression;
  /**
   * @brief Only enabled for @ref basematrix "base matrices" or simple @ref boundarymatrix "boundary matrices", i.e.,
   * when both @ref has_vine_update and @ref can_retrieve_representative_cycles are false.
   * If set to true, the methods @ref Matrix::swap_columns and @ref Matrix::swap_rows are enabled.
   */
  static const bool has_column_and_row_swaps;

  /**
   * @brief If set to true, the underlying container containing the matrix columns is an std::unordered_map. 
   * If set to false, the container is a std::vector. By default, it is recommended to set it to false, but some 
   * methods require it to be true to be enabled: 
   * - @ref Matrix::remove_column(index_type) "remove_column(index)" for @ref basematrix "base matrices",
   * - @ref Matrix::remove_maximal_face(index_type) "remove_maximal_face(index)" for @ref chainmatrix "chain matrices",
   * - @ref Matrix::remove_maximal_face(index_type, const std::vector<index_type>&)
   *  "remove_maximal_face(id_index, const std::vector<id_index>&)" for @ref chainmatrix "chain matrices",
   * - @ref Matrix::remove_last "remove_last()" for @ref chainmatrix "chain matrices" if @ref has_vine_update is true.
   */
  static const bool has_map_column_container;
  /**
   * @brief If set to true, enables the methods @ref Matrix::remove_maximal_face and @ref Matrix::remove_last,
   * except for @ref basematrix "base matrices" when @ref has_column_compression is true.
   */
  static const bool has_removable_columns;

  /**
   * @brief If set to true, enables the method @ref Matrix::get_row.
   */
  static const bool has_row_access;
  /**
   * @brief Only enabled if @ref has_row_access is true, ignored otherwise.
   * If set to true, the underlying container representing a row is an boost::intrusive::list<Cell>.
   * If set to false, the container is a std::set<Cell>. It is usually recommended to set it to true.
   */
  static const bool has_intrusive_rows;
  /**
   * @brief Only enabled if @ref has_row_access is true, ignored otherwise.
   * If set to true, the underlying container containing the rows is an std::map and for
   * @ref chainmatrix "chain matrices", enables the method @ref Matrix::erase_row (always enabled for other
   * matrix types). If set to false, the container is a std::vector.
   */
  static const bool has_removable_rows;

  /**
   * @brief Only used, when at least one of the following is true:
   * @ref has_column_pairings, @ref has_vine_update or @ref can_retrieve_representative_cycles. Is ignored otherwise.
   * If set to true, the matrix is a @ref boundarymatrix "boundary matrix". If set to false, the matrix is a
   * @ref chainmatrix "chain matrix".
   */
  static const bool is_of_boundary_type;

  /**
   * @brief Only enabled for @ref boundarymatrix "boundary" and @ref chainmatrix "chain matrices", i.e., when at least
   * one of the following is true: @ref has_column_pairings, @ref has_vine_update or
   * @ref can_retrieve_representative_cycles. Is ignored otherwise (the notion of dimension makes generally no
   * sense then). If set to true, enables the method @ref Matrix::get_max_dimension. If set to false, the method is
   * disabled except when @ref has_column_pairings is true and @ref has_vine_update and
   * @ref can_retrieve_representative_cycles are both false. In this case, the method is always available.
   */
  static const bool has_matrix_maximal_dimension_access;
  /**
   * @brief If set to true, enables the method @ref Matrix::get_current_barcode. The matrix will then either be a
   * @ref boundarymatrix "boundary matrix" (if @ref is_of_boundary_type is true), or a @ref chainmatrix "chain matrix"
   * (if @ref is_of_boundary_type is false).
   */
  static const bool has_column_pairings;
  /**
   * @brief If set to true, enables the methods @ref Matrix::vine_swap and @ref Matrix::vine_swap_with_z_eq_1_case.
   * The matrix will then either be a @ref boundarymatrix "boundary matrix" (if @ref is_of_boundary_type is true),
   * or a @ref chainmatrix "chain matrix" (if @ref is_of_boundary_type is false).
   */
  static const bool has_vine_update;
  /**
   * @brief If set to true, enables the methods @ref Matrix::update_representative_cycles and
   * @ref Matrix::get_representative_cycles.
   * The matrix will then either be a @ref boundarymatrix "boundary matrix" (if @ref is_of_boundary_type is true),
   * or a @ref chainmatrix "chain matrix" (if @ref is_of_boundary_type is false).
   */
  static const bool can_retrieve_representative_cycles;

  // not implemented yet
  // /**
  //  * @brief Only enabled for boundary and @ref chainmatrix "chain matrices", i.e., when at least one of
  //  * the following is true: @ref has_column_pairings, @ref has_vine_update or
  //  * @ref can_retrieve_representative_cycles.
  //  * Is ignored otherwise
  //  * If set to true, the matrix is decomposed in several submatrices containing each all the
  //  * columns of same dimension.
  //  */
  //  static const bool is_separated_by_dimension;
  //  /**
  //   * @brief If set to true, some methods will use parallel computing.
  //   */
  //  static const bool is_parallelizable;
};

}  // namespace persistence_matrix
}  // namespace Gudhi

