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
 * @file boundary_cell_position_to_id_mapper.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Cell_position_to_ID_mapper class and
 * @ref Gudhi::persistence_matrix::Dummy_pos_mapper structure.
 */

#ifndef PM_ID_POS_MAPPER_H
#define PM_ID_POS_MAPPER_H

#include <unordered_map>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @private
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Cell_position_to_ID_mapper.
 */
struct Dummy_pos_mapper {
  friend void swap([[maybe_unused]] Dummy_pos_mapper& d1, [[maybe_unused]] Dummy_pos_mapper& d2) {}
};

/**
 * @private
 * @ingroup persistence_matrix
 *
 * @brief Map from cell position to cell ID. Only stores a pair if ID != position and has_removable_column is true.
 * 
 * @tparam ID_index @ref IDIdx index type
 * @tparam Pos_index @ref PosIdx index type
 */
template<typename ID_index, typename Pos_index>
struct Cell_position_to_ID_mapper {
  using Index_map = std::unordered_map<Pos_index,ID_index>; //TODO: test other map types

  Index_map map_;

  friend void swap(Cell_position_to_ID_mapper& mapper1, Cell_position_to_ID_mapper& mapper2) {
    mapper1.map_.swap(mapper2.map_);
  }
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_ID_POS_MAPPER_H
