/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_FILTRATION_VALUE_UTILS_H_
#define SIMPLEX_TREE_FILTRATION_VALUE_UTILS_H_

#include <limits>   // std::numeric_limits
#include <cmath>    // std::isnan

#include <gudhi/serialization_utils.h>  // access through this file to serialization methods for arithmetic types

namespace Gudhi {

namespace simplex_tree {

/**
 * @ingroup simplex_tree
 * @brief Returns `Filtration_value(0)` when converted to `Filtration_value`.
 */
inline const struct empty_filtration_value_t {
  template <class T> explicit operator T() const { return T(0); }
} empty_filtration_value;

}  // namespace simplex_tree

/**
 * @ingroup simplex_tree
 * @brief Returns true if and only if the given filtration value is at infinity.
 * This is the overload for when @ref FiltrationValue is an arithmetic type, like double, int etc. It simply
 * tests equality with `std::numeric_limits<FiltrationValue>::infinity()` if defined or with
 * `std::numeric_limits<FiltrationValue>::max()` otherwise. Can therefore be also used with other classes
 * as long as infinity is defined that way.
 */
template <typename Arithmetic_filtration_value>
bool is_positive_infinity(const Arithmetic_filtration_value& f)
{
  if constexpr (std::numeric_limits<Arithmetic_filtration_value>::has_infinity) {
    return f == std::numeric_limits<Arithmetic_filtration_value>::infinity();
  } else {
    return f == std::numeric_limits<Arithmetic_filtration_value>::max();
  }
}

/**
 * @ingroup simplex_tree
 * @brief Given two filtration values at which a simplex exists, stores in the first value the minimal union of births
 * generating a lifetime including those two values.
 * This is the overload for when @ref FiltrationValue is an arithmetic type, like double, int etc.
 * Because the filtration values are totally ordered then, the union is simply the minimum of the two values.
 *
 * NaN values are not supported.
 */
template <typename Arithmetic_filtration_value>
bool unify_lifetimes(Arithmetic_filtration_value& f1, const Arithmetic_filtration_value& f2)
{
  if (f2 < f1){
    f1 = f2;
    return true;
  }
  return false;
}

/**
 * @ingroup simplex_tree
 * @brief Given two filtration values, stores in the first value the lowest common upper bound of the two values.
 * If a filtration value has value `NaN`, it should be considered as the lowest value possible.
 * This is the overload for when @ref FiltrationValue is an arithmetic type, like double, float, int etc.
 * Because the filtration values are totally ordered then, the upper bound is always the maximum of the two values.
 */
template <typename Arithmetic_filtration_value>
bool intersect_lifetimes(Arithmetic_filtration_value& f1, const Arithmetic_filtration_value& f2)
{
  if constexpr (std::numeric_limits<Arithmetic_filtration_value>::has_quiet_NaN) {
    if (std::isnan(f1)) {
      f1 = f2;
      return !std::isnan(f2);
    }

    // Computes the max while handling NaN as lowest value.
    if (!(f1 < f2)) return false;

    f1 = f2;
    return true;
  } else {
    // NaN not possible.
    if (f1 < f2){
      f1 = f2;
      return true;
    }
    return false;
  }
}

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_FILTRATION_VALUE_UTILS_H_
