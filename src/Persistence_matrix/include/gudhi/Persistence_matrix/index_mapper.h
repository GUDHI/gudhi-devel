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
 * @file index_mapper.h
 * @author Hannah Schreiber
 * @brief Contains the Gudhi::persistence_matrix::Index_mapper class and
 * Gudhi::persistence_matrix::Dummy_index_mapper structure.
 */

#ifndef PM_INDEX_MAPPER_H
#define PM_INDEX_MAPPER_H

namespace Gudhi {
namespace persistence_matrix {

// Note: this would be a good candidate for std::optional instead of a mixin as its only role for now
// is to store a map. If this does not change later.

/**
 * @private
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Index_mapper.
 */
struct Dummy_index_mapper {
  friend void swap([[maybe_unused]] Dummy_index_mapper& d1, [[maybe_unused]] Dummy_index_mapper& d2) noexcept {}
};

/**
 * @private
 * @ingroup persistence_matrix
 *
 * @brief Map container. Though for translation between different index types.
 *
 * @tparam Map Map type
 */
template <class Map>
struct Index_mapper {
  using Index_map = Map;

  Index_map map_;

  friend void swap(Index_mapper& mapper1, Index_mapper& mapper2) noexcept { mapper1.map_.swap(mapper2.map_); }
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_INDEX_MAPPER_H
