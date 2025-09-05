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
 * @file chain_pairing.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Chain_pairing class and
 * @ref Gudhi::persistence_matrix::Dummy_chain_pairing structure.
 */

#ifndef PM_CHAIN_PAIRING_H
#define PM_CHAIN_PAIRING_H

namespace Gudhi {
namespace persistence_matrix {

template <typename Master_matrix>
class Chain_barcode_swap;

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Chain_pairing, when the computation of the barcode was not enabled or if the pairing
 * is already managed by the vine update classes.
 */
struct Dummy_chain_pairing {
  friend void swap([[maybe_unused]] Dummy_chain_pairing& d1, [[maybe_unused]] Dummy_chain_pairing& d2) noexcept {}
};

/**
 * @class Chain_pairing chain_pairing.h gudhi/Persistence_matrix/chain_pairing.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the barcode for @ref Chain_matrix if the option was enabled.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Chain_pairing
{
 public:
  using Barcode = typename Master_matrix::Barcode;     /**< Barcode type. */
  using Dimension = typename Master_matrix::Dimension; /**< Dimension value type. */

  /**
   * @brief Default constructor.
   */
  Chain_pairing() = default;

  /**
   * @brief Returns the current barcode which is maintained at any insertion, removal or vine swap.
   *
   * @return Const reference to the barcode.
   */
  const Barcode& get_current_barcode() const;

  /**
   * @brief Swap operator.
   */
  friend void swap(Chain_pairing& pairing1, Chain_pairing& pairing2) noexcept
  {
    pairing1.barcode_.swap(pairing2.barcode_);
    pairing1.indexToBar_.swap(pairing2.indexToBar_);
  }

 protected:
  using Pos_index = typename Master_matrix::Pos_index;
  using Index = typename Master_matrix::Index;

  void _update_barcode(Pos_index birth, Pos_index death);
  void _add_bar(Dimension dim, Pos_index birth);
  void _erase_bar(Pos_index event);
  Pos_index _death(Index index) const;
  Pos_index _birth(Index index) const;
  void _reset();

 private:
  using Dictionary = typename Master_matrix::Bar_dictionary;

  // could also just mark everything as protected as Chain_barcode_swap inherits from Chain_pairing
  // but this way, it marks a better difference between "class using this mixin" with "class extending this mixin"
  friend Chain_barcode_swap<Master_matrix>;

  Barcode barcode_;       /**< Bar container. */
  Dictionary indexToBar_; /**< Map from @ref MatIdx index to bar index. */
};

template <class Master_matrix>
inline const typename Chain_pairing<Master_matrix>::Barcode& Chain_pairing<Master_matrix>::get_current_barcode() const
{
  return barcode_;
}

template <class Master_matrix>
inline void Chain_pairing<Master_matrix>::_update_barcode(Pos_index birth, Pos_index death)
{
  if constexpr (Master_matrix::Option_list::has_removable_columns) {
    auto& barIt = indexToBar_.at(birth);
    barIt->death = death;
    indexToBar_.try_emplace(death, barIt);  // list so iterators are stable
  } else {
    barcode_[indexToBar_[birth]].death = death;
    indexToBar_.push_back(indexToBar_[birth]);
  }
}

template <class Master_matrix>
inline void Chain_pairing<Master_matrix>::_add_bar(Dimension dim, Pos_index birth)
{
  barcode_.emplace_back(birth, Master_matrix::template get_null_value<Pos_index>(), dim);
  if constexpr (Master_matrix::Option_list::has_removable_columns) {
    indexToBar_.try_emplace(birth, --barcode_.end());
  } else {
    indexToBar_.push_back(barcode_.size() - 1);
  }
}

template <class Master_matrix>
inline void Chain_pairing<Master_matrix>::_erase_bar(Pos_index event)
{
  auto it = indexToBar_.find(event);
  typename Barcode::iterator bar = it->second;

  if (bar->death == Master_matrix::template get_null_value<Pos_index>())
    barcode_.erase(bar);
  else
    bar->death = Master_matrix::template get_null_value<Pos_index>();

  indexToBar_.erase(it);
}

template <class Master_matrix>
inline typename Chain_pairing<Master_matrix>::Pos_index Chain_pairing<Master_matrix>::_death(Index index) const
{
  if constexpr (Master_matrix::Option_list::has_removable_columns) {
    return indexToBar_.at(index)->death;
  } else {
    return barcode_.at(indexToBar_.at(index)).death;
  }
}

template <class Master_matrix>
inline typename Chain_pairing<Master_matrix>::Pos_index Chain_pairing<Master_matrix>::_birth(Index index) const
{
  if constexpr (Master_matrix::Option_list::has_removable_columns) {
    return indexToBar_.at(index)->birth;
  } else {
    return barcode_.at(indexToBar_.at(index)).birth;
  }
}

template <class Master_matrix>
inline void Chain_pairing<Master_matrix>::_reset()
{
  barcode_.clear();
  indexToBar_.clear();
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_CHAIN_PAIRING_H
