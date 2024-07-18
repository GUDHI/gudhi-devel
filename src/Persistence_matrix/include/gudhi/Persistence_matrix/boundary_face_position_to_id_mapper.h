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
 * @file base_pairing.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Base_pairing class and
 * @ref Gudhi::persistence_matrix::Dummy_base_pairing structure.
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
 * Inherited instead of @ref Face_position_to_ID_mapper.
 */
struct Dummy_pos_mapper {
  friend void swap([[maybe_unused]] Dummy_pos_mapper& d1, [[maybe_unused]] Dummy_pos_mapper& d2) {}
};

/**
 * @private
 * @ingroup persistence_matrix
 *
 * @brief Map from face position to face ID. Only stores a pair if ID != position and has_removable_column is true.
 * 
 * @tparam id_index @ref IDIdx index type
 * @tparam pos_index @ref PosIdx index type
 */
template<typename id_index, typename pos_index>
struct Face_position_to_ID_mapper {
  using map_type = std::unordered_map<pos_index,id_index>; //TODO: test other map types

  map_type map_;

  friend void swap(Face_position_to_ID_mapper& mapper1, Face_position_to_ID_mapper& mapper2) {
    mapper1.map_.swap(mapper2.map_);
  }
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_ID_POS_MAPPER_H
