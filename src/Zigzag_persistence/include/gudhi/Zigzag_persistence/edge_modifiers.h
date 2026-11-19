/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria and Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file edge_modifiers.h
 * @author Clément Maria, Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::zigzag_persistence::Identity_edge_modifier class and
 * @ref Gudhi::zigzag_persistence::Square_root_edge_modifier class.
 */

#ifndef ZIGZAG_EDGE_MODIFIERS_H_
#define ZIGZAG_EDGE_MODIFIERS_H_

#include <cmath>

#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @class Identity_edge_modifier edge_modifiers.h gudhi/Zigzag_persistence/edge_modifiers.h
 * @brief Identity modifier, i.e., does nothing.
 *
 * @ingroup zigzag_persistence
 */
class Identity_edge_modifier
{
 public:
  /**
   * @brief Returns the given value.
   */
  template <typename Filtration_value>
  static Filtration_value apply_modifier(Filtration_value f)
  {
    return f;
  }

  /**
   * @brief Returns the given value.
   */
  template <typename Filtration_value>
  static Filtration_value apply_inverse_modifier(Filtration_value f)
  {
    return f;
  }

 private:
  /**
   * @brief Default constructor. Should not be called and therefore private.
   */
  Identity_edge_modifier() {}
};

/**
 * @class Square_root_edge_modifier edge_modifiers.h gudhi/Zigzag_persistence/edge_modifiers.h
 * @brief Modifier that square roots the filtration value of an edge.
 *
 * @ingroup zigzag_persistence
 *
 * @details Useful in particular when geometric computations (edge length, etc) are
 * run with squared Euclidean distance for performance.
 *
 * @tparam Filtration_value Filtration value type. Should be compatible with `std::sqrt` and `operator*`.
 */
template <typename Filtration_value>
class Square_root_edge_modifier
{
 public:
  /**
   * @brief Returns the square root of the given value. The parameter it-self is not modified.
   *
   * @param f Value to modify.
   * @return The modified value of @p f.
   */
  static Filtration_value apply_modifier(Filtration_value f) { return std::sqrt(f); }

  /**
   * @brief Returns the square of the given value. The parameter it-self is not modified.
   *
   * @param f Value to modify.
   * @return The modified value of @p f.
   */
  static Filtration_value apply_inverse_modifier(Filtration_value f) { return f * f; }

 private:
  /**
   * @brief Default constructor. Should not be called and therefore private.
   */
  Square_root_edge_modifier() {}
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // ZIGZAG_EDGE_MODIFIERS_H_
