/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ZZ_OPTIONS_TYPE_H_
#define CONCEPT_ZZ_OPTIONS_TYPE_H_

/** @file ZigzagPersistenceOptions.h
 * @brief Contains @ref Gudhi::zigzag_persistence::ZigzagPersistenceOptions concept.
 */

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief List of options used for the matrix maintained for the zigzag persistence computation.
 */
struct ZigzagPersistenceOptions {
  /**
   * @brief Type for the coefficient field type. Has to support \f$Z_2\f$.
   */
  typename field_coeff_type;

  /**
   * @brief Has to be set to true. Indicates that the computation will be made with \f$Z_2\f$ coefficients.
   */
  static const bool is_z2 = true;
  /**
   * @brief Type of the columns in the matrix. 
   * The available column types are given by @ref Gudhi::persistence_matrix::Column_types. 
   * The column type has to support row access.
   */
  static const Column_types column_type;

  /**
   * @brief Has to be set to true. Indicates that the rows should be directly accessible in the matrix.
   */
  static const bool has_row_access = true;
  /**
   * @brief Set to true, if the rows should be intrusive lists or to false if they should be sets. True is recommended.
   * Note that intrusive rows are not compatible with certain column types.
   */
  static const bool has_intrusive_rows;
  /**
   * @brief Has to set to true. Indicates that the rows of the matrix can be removed.
   */
  static const bool has_removable_rows = true;
  /**
   * @brief Has to be set to false. Indicates that the matrix should not store birth/death pairs of its columns.
   */
  static const bool has_column_pairings = false;
  /**
   * @brief Has to be set to true. Enables maintaining the matrix while switching columns.
   */
  static const bool has_vine_update = true;
  /**
   * @brief If set to true, the matrix can retrieve the representative cycles for the cycle classes. 
   * This option is useless for zigzag computation and therefore it is recommended to set it to false.
   */
  static const bool can_retrieve_representative_cycles;
  /**
   * @brief This value has to be defined but will be ignored.
   */
  static const bool has_column_compression;
  /**
   * @brief Has to be set to false. 
   * Indicates that the matrix should represent the base of the chain complex and not of the boundary group.
   */
  static const bool is_of_boundary_type = false;
  /**
   * @brief Has to be set to true. Indicates that the columns of the matrix can be removed.
   */
  static const bool has_removable_columns = true;
  /**
   * @brief Has to be set to false.
   * Indicates that the access to the columns will be done through simplex IDs instead of column positions.
   */
  static const bool is_indexed_by_position = false;
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // CONCEPT_ZZ_OPTIONS_TYPE_H_
