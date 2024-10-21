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
 * @file ru_pairing.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::RU_pairing class and
 * @ref Gudhi::persistence_matrix::Dummy_ru_pairing structure.
 */

#ifndef PM_RU_PAIRING_H
#define PM_RU_PAIRING_H

#include <unordered_map>

#include "boundary_cell_position_to_id_mapper.h"

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref RU_pairing, when the computation of the barcode was not enabled or if the pairing
 * is already managed by the vine update classes.
 */
struct Dummy_ru_pairing
{
  friend void swap([[maybe_unused]] Dummy_ru_pairing& d1, [[maybe_unused]] Dummy_ru_pairing& d2) {}
};

/**
 * @class RU_pairing ru_pairing.h gudhi/Persistence_matrix/ru_pairing.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the barcode for @ref RU_matrix if the option was enabled.
 * 
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class RU_pairing : public std::conditional<
                       Master_matrix::Option_list::has_removable_columns,
                       Cell_position_to_ID_mapper<typename Master_matrix::ID_index, typename Master_matrix::Pos_index>,
                       Dummy_pos_mapper
                    >::type
{
 protected:
  using Pos_index = typename Master_matrix::Pos_index;
  using ID_index = typename Master_matrix::ID_index;
  //PIDM = Position to ID Map
  using PIDM = typename std::conditional<Master_matrix::Option_list::has_removable_columns,
                                          Cell_position_to_ID_mapper<ID_index, Pos_index>,
                                          Dummy_pos_mapper
                                        >::type;

 public:
  using Barcode = typename Master_matrix::Barcode;  /**< Barcode type. */

  /**
   * @brief Default constructor.
   */
  RU_pairing() : PIDM() {}

  /**
   * @brief Returns the current barcode which is maintained at any insertion, removal or vine swap.
   * 
   * @return Const reference to the barcode.
   */
  const Barcode& get_current_barcode() const { return barcode_; }

  /**
   * @brief Swap operator.
   */
  friend void swap(RU_pairing& pairing1, RU_pairing& pairing2) {
    swap(static_cast<PIDM&>(pairing1), static_cast<PIDM&>(pairing2));
    pairing1.barcode_.swap(pairing2.barcode_);
    pairing1.indexToBar_.swap(pairing2.indexToBar_);
    pairing1.idToPosition_.swap(pairing2.idToPosition_);
  }

 protected:
  using Dimension = typename Master_matrix::Dimension;
  using Dictionary = typename Master_matrix::Bar_dictionary;

  Barcode barcode_;       /**< Bar container. */
  Dictionary indexToBar_; /**< Map from @ref MatIdx index to bar index. */
  /**
   * @brief Map from cell ID to cell position. Only stores a pair if ID != position.
   */
  std::unordered_map<ID_index, Pos_index> idToPosition_;  //TODO: test other map types

  void _update_barcode(ID_index birthPivot, Pos_index death) {
    auto it = idToPosition_.find(birthPivot);
    Pos_index pivotBirth = it == idToPosition_.end() ? birthPivot : it->second;
    if constexpr (Master_matrix::hasFixedBarcode || !Master_matrix::Option_list::has_removable_columns) {
      barcode_[indexToBar_[pivotBirth]].death = death;
      indexToBar_.push_back(indexToBar_[pivotBirth]);
    } else {
      auto& barIt = indexToBar_.at(pivotBirth);
      barIt->death = death;
      indexToBar_.try_emplace(death, barIt);  // list so iterators are stable
    }
  }

  void _add_bar(Dimension dim, Pos_index birth) {
    barcode_.emplace_back(birth, -1, dim);
    if constexpr (Master_matrix::hasFixedBarcode || !Master_matrix::Option_list::has_removable_columns) {
      indexToBar_.push_back(barcode_.size() - 1);
    } else {
      indexToBar_.try_emplace(birth, --barcode_.end());
    }
  }

  void _remove_last(Pos_index eventIndex) {
    static_assert(Master_matrix::Option_list::has_removable_columns, "_remove_last not available.");

    if constexpr (Master_matrix::hasFixedBarcode) {
      auto& bar = barcode_[indexToBar_[eventIndex]];
      if (bar.death == static_cast<Pos_index>(-1)) {  // birth
        barcode_.pop_back();  // sorted by birth and eventIndex has to be the highest one
      } else {                // death
        bar.death = -1;
      };
      indexToBar_.pop_back();
    } else {  // birth order eventually shuffled by vine updates. No sort possible to keep the matchings.
      auto it = indexToBar_.find(eventIndex);
      typename Barcode::iterator bar = it->second;

      if (bar->death == static_cast<Pos_index>(-1))
        barcode_.erase(bar);
      else
        bar->death = -1;

      indexToBar_.erase(it);
    }

    auto it = PIDM::map_.find(eventIndex);
    if (it != PIDM::map_.end()){
      idToPosition_.erase(it->second);
      PIDM::map_.erase(it);
    }
  }
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_RU_PAIRING_H
