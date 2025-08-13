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
 * @file chain_vine_swap.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Chain_barcode_swap and
 * @ref Gudhi::persistence_matrix::Chain_vine_swap classes, as well as the
 * @ref Gudhi::persistence_matrix::Dummy_chain_vine_swap and
 * @ref Gudhi::persistence_matrix::Dummy_chain_vine_pairing structures.
 */

#ifndef PM_CHAIN_VINE_SWAP_H
#define PM_CHAIN_VINE_SWAP_H

#include <cassert>
#include <utility>     //std::swap & std::move
#include <functional>  //std::function
#include <stdexcept>   //std::invalid_argument

#include <gudhi/Persistence_matrix/chain_pairing.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Default death comparator. Simply assumes that two positive paired columns are never swapped. Which is true
 * for the use case in zigzag persistence for example.
 *
 * @param columnIndex1 First column index.
 * @param columnIndex2 Second column index.
 * @return false
 */
constexpr bool _no_G_death_comparator([[maybe_unused]] unsigned int columnIndex1,
                                      [[maybe_unused]] unsigned int columnIndex2)
{
  return false;
}

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Chain_vine_swap, when vine swaps are not enabled.
 */
struct Dummy_chain_vine_swap {
  using CP = void;

  friend void swap([[maybe_unused]] Dummy_chain_vine_swap& d1, [[maybe_unused]] Dummy_chain_vine_swap& d2) noexcept {}

  Dummy_chain_vine_swap() = default;

  template <typename BirthComparatorFunction, typename DeathComparatorFunction>
  Dummy_chain_vine_swap([[maybe_unused]] const BirthComparatorFunction& birthComparator,
                        [[maybe_unused]] const DeathComparatorFunction& deathComparator)
  {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Chain_barcode_swap, when the barcode is not stored.
 */
struct Dummy_chain_vine_pairing {
  friend void swap([[maybe_unused]] Dummy_chain_vine_pairing& d1,
                   [[maybe_unused]] Dummy_chain_vine_pairing& d2) noexcept
  {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Class managing the barcode for @ref Chain_vine_swap.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <typename Master_matrix>
class Chain_barcode_swap : public Chain_pairing<Master_matrix>
{
 public:
  using ID_index = typename Master_matrix::ID_index;   /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index; /**< @ref PosIdx index type. */
  // CP = Chain Pairing
  using CP = Chain_pairing<Master_matrix>;

  /**
   * @brief Default constructor.
   */
  Chain_barcode_swap() = default;
  /**
   * @brief Copy constructor.
   *
   * @param toCopy Matrix to copy.
   */
  Chain_barcode_swap(const Chain_barcode_swap& toCopy) : CP(static_cast<const CP&>(toCopy)) {};
  /**
   * @brief Move constructor.
   *
   * @param other Matrix to move.
   */
  Chain_barcode_swap(Chain_barcode_swap&& other) noexcept : CP(std::move(static_cast<CP&>(other))) {};

  ~Chain_barcode_swap() = default;

  Chain_barcode_swap& operator=(Chain_barcode_swap other)
  {
    CP::operator=(other);
    return *this;
  }

  Chain_barcode_swap& operator=(Chain_barcode_swap&& other) noexcept
  {
    CP::operator=(std::move(other));
    return *this;
  }

  friend void swap(Chain_barcode_swap& swap1, Chain_barcode_swap& swap2) noexcept
  {
    swap(static_cast<CP&>(swap1), static_cast<CP&>(swap2));
  }

 protected:
  bool _is_negative_in_bar(ID_index pivot) const
  {
    Pos_index pos = _get_pivot_position(pivot);
    return _death_val(pivot) == pos;
  }

  void _positive_transpose_barcode(ID_index pivot1, ID_index pivot2)
  {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);

    _birth(pos1) = pos2;
    _birth(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  void _negative_transpose_barcode(ID_index pivot1, ID_index pivot2)
  {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);

    _death(pos1) = pos2;
    _death(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  void _positive_negative_transpose_barcode(ID_index pivot1, ID_index pivot2)
  {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);

    _birth(pos1) = pos2;
    _death(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  void _negative_positive_transpose_barcode(ID_index pivot1, ID_index pivot2)
  {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);

    _death(pos1) = pos2;
    _birth(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  bool _are_adjacent(ID_index pivot1, ID_index pivot2) const
  {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);
    return pos1 < pos2 ? (pos2 - pos1) == 1 : (pos1 - pos2) == 1;
  }

  Pos_index _death_val(ID_index pivot) const { return CP::_death(_get_pivot_position(pivot)); }

  Pos_index _birth_val(ID_index pivot) const { return CP::_birth(_get_pivot_position(pivot)); }

  void _reset() { CP::_reset(); }

 private:
  using Master_chain_matrix = typename Master_matrix::Master_chain_matrix;

  Pos_index _get_pivot_position(ID_index pivot) const
  {
    const auto& map = static_cast<const Master_chain_matrix*>(this)->map_;
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      // TODO: quite often called, make public and pass position instead of pivot to avoid find() every time?
      return map.at(pivot);
    } else {
      return map[pivot];
    }
  }

  Pos_index& _death(Pos_index index)
  {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return CP::indexToBar_.at(index)->death;
    } else {
      return CP::barcode_.at(CP::indexToBar_.at(index)).death;
    }
  }

  Pos_index& _birth(Pos_index index)
  {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return CP::indexToBar_.at(index)->birth;
    } else {
      return CP::barcode_.at(CP::indexToBar_.at(index)).birth;
    }
  }
};

/**
 * @class Chain_vine_swap chain_vine_swap.h gudhi/Persistence_matrix/chain_vine_swap.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the vine swaps for @ref Chain_matrix.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Chain_vine_swap
{
 public:
  using Index = typename Master_matrix::Index;                       /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;                 /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index;               /**< @ref PosIdx index type. */
  using Column_container = typename Master_matrix::Column_container; /**< Column container type. */
  using Column = typename Master_matrix::Column;                     /**< Column type. */
  using EventCompFuncPointer = bool (*)(Pos_index, Pos_index);       /**< Pointer type for birth/death comparators. */

  /**
   * @brief Default constructor. Only available if @ref PersistenceMatrixOptions::has_column_pairings is true.
   */
  Chain_vine_swap();
  /**
   * @brief Constructor storing the given comparators.
   *
   * @param birthComparator Method taking two @ref PosIdx indices as input and returning true if and only if
   * the birth associated to the first position is strictly less than birth associated to
   * the second one with respect to some self defined order. It is used while swapping two unpaired or
   * two negative columns.
   * @param deathComparator Method taking two @ref PosIdx indices as input and returning true if and only if
   * the death associated to the first position is strictly less than death associated to
   * the second one with respect to some self defined order. It is used while swapping two positive but paired
   * columns. Default value: @ref _no_G_death_comparator.
   */
  Chain_vine_swap(std::function<bool(Pos_index, Pos_index)> birthComparator,
                  std::function<bool(Pos_index, Pos_index)> deathComparator = _no_G_death_comparator);

  /**
   * @brief Does the same than @ref vine_swap, but assumes that the swap is non trivial and
   * therefore skips a part of the case study.
   *
   * @param columnIndex1 @ref MatIdx index of the first cell.
   * @param columnIndex2 @ref MatIdx index of the second cell.
   * @return Let \f$ pos1 \f$ be the @ref PosIdx index of @p columnIndex1 and \f$ pos2 \f$ be the @ref PosIdx index of
   * @p columnIndex2. The method returns the @ref MatIdx of the column which has now, after the swap, the @ref PosIdx
   * \f$ max(pos1, pos2) \f$.
   */
  Index vine_swap_with_z_eq_1_case(Index columnIndex1, Index columnIndex2);
  /**
   * @brief Does a vine swap between two cells which are consecutive in the filtration. Roughly, if \f$ F \f$ is
   * the current filtration represented by the matrix, the method modifies the matrix such that the new state
   * corresponds to a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but with the two given cells
   * at swapped positions. Of course, the two cells should not have a face/coface relation which each other ;
   * \f$ F' \f$ has to be a valid filtration.
   * See @cite vineyards for more information about vine and vineyards.
   *
   * @param columnIndex1 @ref MatIdx index of the first cell.
   * @param columnIndex2 @ref MatIdx index of the second cell. It is assumed that the @ref PosIdx of both only differs
   * by one if the barcode is maintained.
   * @return Let \f$ pos1 \f$ be the @ref PosIdx index of @p columnIndex1 and \f$ pos2 \f$ be the @ref PosIdx index of
   * @p columnIndex2. The method returns the @ref MatIdx of the column which has now, after the swap, the @ref PosIdx
   * \f$ max(pos1, pos2) \f$.
   */
  Index vine_swap(Index columnIndex1, Index columnIndex2);

  /**
   * @brief Swap operator.
   */
  friend void swap(Chain_vine_swap& swap1, Chain_vine_swap& swap2) noexcept
  {
    std::swap(swap1.birthComp_, swap2.birthComp_);
    std::swap(swap1.deathComp_, swap2.deathComp_);
  }

 private:
  using Master_chain_matrix = typename Master_matrix::Master_chain_matrix;

  std::function<bool(Pos_index, Pos_index)> birthComp_; /**< for F x F & H x H. */
  std::function<bool(Pos_index, Pos_index)> deathComp_; /**< for G x G. */

  bool _is_negative_in_pair(Index columnIndex);
  Index _positive_vine_swap(Index columnIndex1, Index columnIndex2);
  Index _positive_negative_vine_swap(Index columnIndex1, Index columnIndex2);
  Index _negative_positive_vine_swap(Index columnIndex1, Index columnIndex2);
  Index _negative_vine_swap(Index columnIndex1, Index columnIndex2);
  void _swap_positions(ID_index pivot1, ID_index pivot2);

  constexpr Master_chain_matrix* _matrix() { return static_cast<Master_chain_matrix*>(this); }

  constexpr const Master_chain_matrix* _matrix() const { return static_cast<const Master_chain_matrix*>(this); }
};

template <class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap() : birthComp_(), deathComp_()
{
  static_assert(Master_matrix::Option_list::has_column_pairings,
                "If barcode is not stored, at least a birth comparator has to be specified.");
}

template <class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(std::function<bool(Pos_index, Pos_index)> birthComparator,
                                                       std::function<bool(Pos_index, Pos_index)> deathComparator)
    : birthComp_(std::move(birthComparator)), deathComp_(std::move(deathComparator))
{}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(
    Index columnIndex1,
    Index columnIndex2)
{
  const bool col1IsNeg = _is_negative_in_pair(columnIndex1);
  const bool col2IsNeg = _is_negative_in_pair(columnIndex2);

  if (col1IsNeg && col2IsNeg) return _negative_vine_swap(columnIndex1, columnIndex2);

  if (col1IsNeg) return _negative_positive_vine_swap(columnIndex1, columnIndex2);

  if (col2IsNeg) return _positive_negative_vine_swap(columnIndex1, columnIndex2);

  return _positive_vine_swap(columnIndex1, columnIndex2);
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::vine_swap(Index columnIndex1,
                                                                                                Index columnIndex2)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    GUDHI_CHECK(_matrix()->_are_adjacent(_matrix()->get_pivot(columnIndex1), _matrix()->get_pivot(columnIndex2)),
                std::invalid_argument(
                    "Chain_vine_swap::vine_swap - Columns to be swapped need to be adjacent in the 'real' matrix."));
  }

  const bool col1IsNeg = _is_negative_in_pair(columnIndex1);
  const bool col2IsNeg = _is_negative_in_pair(columnIndex2);

  if (col1IsNeg && col2IsNeg) {
    if (_matrix()->is_zero_entry(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
      ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
      ID_index pivot2 = _matrix()->get_pivot(columnIndex2);
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        _matrix()->_negative_transpose_barcode(pivot1, pivot2);
      }
      _swap_positions(pivot1, pivot2);
      return columnIndex1;
    }
    return _negative_vine_swap(columnIndex1, columnIndex2);
  }

  if (col1IsNeg) {
    if (_matrix()->is_zero_entry(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
      ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
      ID_index pivot2 = _matrix()->get_pivot(columnIndex2);
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        _matrix()->_negative_positive_transpose_barcode(pivot1, pivot2);
      }
      _swap_positions(pivot1, pivot2);
      return columnIndex1;
    }
    return _negative_positive_vine_swap(columnIndex1, columnIndex2);
  }

  if (col2IsNeg) {
    if (_matrix()->is_zero_entry(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
      ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
      ID_index pivot2 = _matrix()->get_pivot(columnIndex2);
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        _matrix()->_positive_negative_transpose_barcode(pivot1, pivot2);
      }
      _swap_positions(pivot1, pivot2);
      return columnIndex1;
    }
    return _positive_negative_vine_swap(columnIndex1, columnIndex2);
  }

  if (_matrix()->is_zero_entry(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
    ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
    ID_index pivot2 = _matrix()->get_pivot(columnIndex2);
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      _matrix()->_positive_transpose_barcode(pivot1, pivot2);
    }
    _swap_positions(pivot1, pivot2);
    return columnIndex1;
  }
  return _positive_vine_swap(columnIndex1, columnIndex2);
}

template <class Master_matrix>
inline bool Chain_vine_swap<Master_matrix>::_is_negative_in_pair(Index columnIndex)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    return _matrix()->_is_negative_in_bar(_matrix()->get_pivot(columnIndex));
  } else {
    auto& col = _matrix()->get_column(columnIndex);
    if (!col.is_paired()) return false;
    return col.get_pivot() > _matrix()->get_pivot(col.get_paired_chain_index());
  }
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::_positive_vine_swap(
    Index columnIndex1,
    Index columnIndex2)
{
  auto& col1 = _matrix()->get_column(columnIndex1);
  auto& col2 = _matrix()->get_column(columnIndex2);

  _swap_positions(col1.get_pivot(), col2.get_pivot());
  // TODO: factorize the cases
  // But for debug it is much more easier to understand what is happening when split like this
  if (!col1.is_paired()) {  // F x *
    bool hasSmallerBirth;
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      // this order because position were swapped with swap_positions
      hasSmallerBirth = (_matrix()->_birth_val(col2.get_pivot()) < _matrix()->_birth_val(col1.get_pivot()));
    } else {
      hasSmallerBirth = birthComp_(columnIndex1, columnIndex2);
    }

    if (!col2.is_paired() && hasSmallerBirth) {
      _matrix()->add_to(columnIndex1, columnIndex2);
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        _matrix()->_positive_transpose_barcode(col1.get_pivot(), col2.get_pivot());
      }
      return columnIndex1;
    }
    _matrix()->add_to(columnIndex2, columnIndex1);

    return columnIndex2;
  }

  if (!col2.is_paired()) {  // G x F
    static_cast<Master_chain_matrix*>(this)->add_to(columnIndex1, columnIndex2);
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      _matrix()->_positive_transpose_barcode(col1.get_pivot(), col2.get_pivot());
    }
    return columnIndex1;
  }

  bool hasSmallerDeath;
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    // this order because position were swapped with swap_positions
    hasSmallerDeath = (_matrix()->_death_val(col2.get_pivot()) < _matrix()->_death_val(col1.get_pivot()));
  } else {
    hasSmallerDeath = deathComp_(columnIndex1, columnIndex2);
  }

  // G x G
  if (hasSmallerDeath) {
    _matrix()->add_to(col1.get_paired_chain_index(), col2.get_paired_chain_index());
    _matrix()->add_to(columnIndex1, columnIndex2);
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      _matrix()->_positive_transpose_barcode(col1.get_pivot(), col2.get_pivot());
    }
    return columnIndex1;
  }

  _matrix()->add_to(col2.get_paired_chain_index(), col1.get_paired_chain_index());
  _matrix()->add_to(columnIndex2, columnIndex1);

  return columnIndex2;
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::_positive_negative_vine_swap(
    Index columnIndex1,
    Index columnIndex2)
{
  _matrix()->add_to(columnIndex1, columnIndex2);

  ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
  ID_index pivot2 = _matrix()->get_pivot(columnIndex2);
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _matrix()->_positive_negative_transpose_barcode(pivot1, pivot2);
  }
  _swap_positions(pivot1, pivot2);

  return columnIndex1;
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::_negative_positive_vine_swap(
    Index columnIndex1,
    Index columnIndex2)
{
  _matrix()->add_to(columnIndex2, columnIndex1);

  _swap_positions(_matrix()->get_pivot(columnIndex1), _matrix()->get_pivot(columnIndex2));

  return columnIndex2;
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::_negative_vine_swap(
    Index columnIndex1,
    Index columnIndex2)
{
  auto& col1 = _matrix()->get_column(columnIndex1);
  auto& col2 = _matrix()->get_column(columnIndex2);

  Index pairedIndex1 = col1.get_paired_chain_index();
  Index pairedIndex2 = col2.get_paired_chain_index();

  bool hasSmallerBirth;
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    hasSmallerBirth = (_matrix()->_birth_val(col1.get_pivot()) < _matrix()->_birth_val(col2.get_pivot()));
  } else {
    hasSmallerBirth = birthComp_(pairedIndex1, pairedIndex2);
  }

  _swap_positions(col1.get_pivot(), col2.get_pivot());

  if (hasSmallerBirth) {
    _matrix()->add_to(pairedIndex1, pairedIndex2);
    _matrix()->add_to(columnIndex1, columnIndex2);

    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      _matrix()->_negative_transpose_barcode(col1.get_pivot(), col2.get_pivot());
    }

    return columnIndex1;
  }

  _matrix()->add_to(pairedIndex2, pairedIndex1);
  _matrix()->add_to(columnIndex2, columnIndex1);

  return columnIndex2;
}

template <class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_swap_positions(ID_index pivot1, ID_index pivot2)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings ||
                Master_matrix::Option_list::can_retrieve_representative_cycles) {
    auto& map = _matrix()->map_;
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      std::swap(map.at(pivot1), map.at(pivot2));
    } else {
      std::swap(map[pivot1], map[pivot2]);
    }
  }
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_CHAIN_VINE_SWAP_H
