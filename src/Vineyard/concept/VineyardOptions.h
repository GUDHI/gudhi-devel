/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file VineyardOptions.h
 * @author Hannah Schreiber
 * @brief Contains the concept for the vineyard options.
 */

/// Gudhi namespace.
namespace Gudhi {
/// Vineyard namespace.
namespace vineyard {

/**
 * @ingroup vineyard
 *
 * @brief Concept of the template parameter for the class @ref Vineyard_base and @ref Vineyard_builder.
 *
 * An implementation of this concept is @ref Default_vineyard_options.
 * If you want to provide your own, it is recommended that you derive from it and override some parts instead of
 * writing a class from scratch (unless nothing remains common).
 */
struct VineyardOptions {
  /**
   * @brief Type for the dimension. Has to be an integer type.
   * If unsigned, the maximal value of the type should not be attained during a run.
   */
  using Dimension = unspecified;
  /**
   * @brief Type for the different indexation types and should be able to contain the maximal number of cells during
   * a run. Has to be an integer type.
   * If unsigned, the maximal value of the type should not be attained during a run.
   */
  using Index = unspecified;

  /**
   * @brief Indicates the underlying matrix type: either the RU decomposition or a chain complex base. Depending
   * on the data properties one can be faster than the other or vice versa.
   */
  static constexpr bool is_RU;
  /**
   * @brief Specifies the desired column type for the underlying matrix. All possible column types are described in
   * @ref Gudhi::persistence_matrix::Column_types.
   */
  static const Column_types column_type;
};

}  // namespace vineyard
}  // namespace Gudhi
