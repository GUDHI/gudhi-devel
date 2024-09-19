/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ZZ_OPTIONS_TYPE_H_
#define CONCEPT_ZZ_OPTIONS_TYPE_H_

/** @file ZigzagOptions.h
 * @brief Contains @ref Gudhi::zigzag_persistence::ZigzagOptions and
 * @ref Gudhi::zigzag_persistence::FilteredZigzagOptions concept.
 */

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @ingroup zigzag_persistence
 *
 * @brief List of options used for the filtered zigzag persistence computation.
 */
struct FilteredZigzagOptions {
  /**
   * @brief Numerical type for the face IDs used internally and other indexations. It must be signed.
   */
  using Internal_key = unspecified;

  /**
   * @brief Type for the face IDs used at insertion and in the boundaries given as argument.
   * Has to be usable as key in a hashtable, so "hashable" and comparable.
   */
  using Face_key = unspecified;

  /**
   * @brief Type for filtration values.
   */
  using Filtration_value = unspecified;

  /**
   * @brief Type for the dimension values.
   */
  using Dimension = unspecified;

  /**
   * @brief Column type used by the internal matrix.
   */
  static const Gudhi::persistence_matrix::Column_types column_type;
};

/**
 * @ingroup zigzag_persistence
 *
 * @brief List of options used for the zigzag persistence computation.
 */
struct ZigzagOptions {
  /**
   * @brief Numerical type for the face IDs used internally and other indexations. It must be signed.
   */
  using Internal_key = unspecified;

  /**
   * @brief Type for the dimension values.
   */
  using Dimension = unspecified;

  /**
   * @brief Column type used by the internal matrix.
   */
  static const Gudhi::persistence_matrix::Column_types column_type;
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // CONCEPT_ZZ_OPTIONS_TYPE_H_
