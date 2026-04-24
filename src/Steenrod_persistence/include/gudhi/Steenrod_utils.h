/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Anibal M. Medina-Mardones
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef STEENROD_UTILS_H_
#define STEENROD_UTILS_H_

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include <gudhi/Steenrod_types.h>

namespace Gudhi {

namespace steenrod_persistence {

/** \brief Return the symmetric difference of two sorted columns.
 *
 * \ingroup steenrod_persistence
 *
 * This is the GF(2) column-addition used throughout the reduction. It is the
 * innermost operation in both twist reduction and the augmented reduction for
 * Steenrod barcodes.
 *
 * @param[in] a First column, sorted in ascending order.
 * @param[in] b Second column, sorted in ascending order.
 * @return Sorted symmetric difference of ``a`` and ``b``.
 */
inline Column symm_diff(const Column& a, const Column& b) {
  Column result;
  result.reserve(a.size() + b.size());  // worst-case upper bound
  std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result));
  return result;
}

/** \brief In-place variant of \ref symm_diff.
 *
 * \ingroup steenrod_persistence
 *
 * Computes ``a = symm_diff(a, b)`` without requiring the caller to keep the
 * previous contents of ``a``.
 *
 * @param[in,out] a Left operand, overwritten with the symmetric difference.
 * @param[in]     b Right operand, must remain valid.
 */
inline void symm_diff_inplace(Column& a, const Column& b) {
  Column result;
  result.reserve(a.size() + b.size());
  std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result));
  a = std::move(result);
}

/** \brief Return the pivot of a column.
 *
 * \ingroup steenrod_persistence
 *
 * Columns are stored in ascending order, so the pivot is the first element.
 * Returns ``-1`` for an empty column.
 */
inline Index pivot(const Column& col) {
  return col.empty() ? Index{-1} : col.front();
}

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // STEENROD_UTILS_H_
