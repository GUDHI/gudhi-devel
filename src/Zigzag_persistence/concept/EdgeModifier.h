/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ZZ_EDGE_MODIFIER_H_
#define CONCEPT_ZZ_EDGE_MODIFIER_H_

/** @file EdgeFiltrationTransformer.h
 * @brief Contains @ref Gudhi::zigzag_persistence::EdgeFiltrationTransformer concept.
 */

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief Methods whose purposes are to modify the filtration value of a given edge following a rule.
 * The concept is for example realized by @ref Identity_edge_modifier or @ref Square_root_edge_modifier "".
 */
template <typename Filtration_value>
class EdgeFiltrationTransformer {
 public:
  /**
   * @brief Applies the modifier to the given value and returns it.
   * 
   * @param f Value to modify.
   * @return The modified value of @p f.
   */
  static Filtration_value apply_modifier(Filtration_value f);

  /**
   * @brief Applies the inverse modifier to the given value and returns it.
   * So, apply_inverse_modifier(apply_modifier(f)) == f (modulo some possible precision errors.).
   * 
   * @param f Value to modify.
   * @return The modified value of @p f.
   */
  static Filtration_value apply_inverse_modifier(Filtration_value f);
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // CONCEPT_ZZ_EDGE_MODIFIER_H_
