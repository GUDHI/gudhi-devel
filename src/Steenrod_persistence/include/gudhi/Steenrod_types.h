/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Anibal M. Medina-Mardones
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef STEENROD_TYPES_H_
#define STEENROD_TYPES_H_

#include <array>
#include <cstdint>
#include <vector>

namespace Gudhi {

namespace steenrod_persistence {

/** \brief Index type used throughout the module.
 *
 * ``int32_t`` suffices for any realistic filtered complex: simplex indices
 * never exceed 2B, and the 32-bit width improves cache utilisation in the
 * inner reduction loops compared to a 64-bit alternative.
 *
 * \ingroup steenrod_persistence
 */
using Index = int32_t;

/** \brief A column of the boundary/coboundary matrix.
 *
 * Stored as a sorted (ascending) list of row indices. Sorted order is an
 * invariant maintained by every operation acting on a ``Column``.
 *
 * \ingroup steenrod_persistence
 */
using Column = std::vector<Index>;

/** \brief A simplex: a sorted list of vertex indices.
 *
 * \ingroup steenrod_persistence
 */
using Simplex = std::vector<Index>;

/** \brief One persistence bar.
 *
 * Layout is ``{death_index, birth_index}`` (cohomology convention).
 * A ``death_index`` equal to ``-1`` marks an essential (infinite) bar.
 *
 * \ingroup steenrod_persistence
 */
using Bar = std::array<Index, 2>;

/** \brief All bars in one homological dimension. \ingroup steenrod_persistence */
using Barcode = std::vector<Bar>;

/** \brief Full barcode: one ``Barcode`` per dimension. \ingroup steenrod_persistence */
using BarcodeByDim = std::vector<Barcode>;

/** \brief A cohomology representative: a sorted list of *relative* (within-dim) indices.
 * \ingroup steenrod_persistence
 */
using CohoRep = std::vector<Index>;

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // STEENROD_TYPES_H_
