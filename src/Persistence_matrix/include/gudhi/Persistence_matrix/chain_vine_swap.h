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

#include <utility>  //std::swap & std::move
#include <cassert>
#include <functional> //std::function
#include <stdexcept>  //std::invalid_argument

#include "chain_pairing.h"

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
struct Dummy_chain_vine_swap
{
  friend void swap([[maybe_unused]] Dummy_chain_vine_swap& d1, [[maybe_unused]] Dummy_chain_vine_swap& d2) {}

  Dummy_chain_vine_swap() {}
  template <typename BirthComparatorFunction, typename DeathComparatorFunction>
  Dummy_chain_vine_swap([[maybe_unused]] const BirthComparatorFunction& birthComparator,
                        [[maybe_unused]] const DeathComparatorFunction& deathComparator) {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Chain_barcode_swap, when the barcode is not stored.
 */
struct Dummy_chain_vine_pairing
{
  friend void swap([[maybe_unused]] Dummy_chain_vine_pairing& d1, [[maybe_unused]] Dummy_chain_vine_pairing& d2) {}
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
  using ID_index = typename Master_matrix::ID_index;    /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index;  /**< @ref PosIdx index type. */
  //CP = Chain Pairing
  using CP = Chain_pairing<Master_matrix>;

  /**
   * @brief Default constructor.
   */
  Chain_barcode_swap(){};
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Matrix to copy.
   */
  Chain_barcode_swap(const Chain_barcode_swap& toCopy)
      : CP(static_cast<const CP&>(toCopy)), pivotToPosition_(toCopy.pivotToPosition_){};
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Chain_barcode_swap(Chain_barcode_swap&& other)
      : CP(std::move(static_cast<CP&>(other))), pivotToPosition_(std::move(other.pivotToPosition_)){};

 protected:
  using Dictionary = typename Master_matrix::template Dictionary<Pos_index>;

  Dictionary pivotToPosition_;  // necessary to keep track of the barcode changes

  void swap_positions(ID_index pivot1, ID_index pivot2) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      std::swap(pivotToPosition_.at(pivot1), pivotToPosition_.at(pivot2));
    } else {
      std::swap(pivotToPosition_[pivot1], pivotToPosition_[pivot2]);
    }
  }

  bool is_negative_in_pair(ID_index pivot) const {
    Pos_index pos = _get_pivot_position(pivot);
    return death(pivot) == pos;
  }

  void positive_transpose(ID_index pivot1, ID_index pivot2) {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);

    _birth(pos1) = pos2;
    _birth(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  void negative_transpose(ID_index pivot1, ID_index pivot2) {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);

    _death(pos1) = pos2;
    _death(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  void positive_negative_transpose(ID_index pivot1, ID_index pivot2) {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);

    _birth(pos1) = pos2;
    _death(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  void negative_positive_transpose(ID_index pivot1, ID_index pivot2) {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);

    _death(pos1) = pos2;
    _birth(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  bool are_adjacent(ID_index pivot1, ID_index pivot2) const {
    Pos_index pos1 = _get_pivot_position(pivot1);
    Pos_index pos2 = _get_pivot_position(pivot2);
    return pos1 < pos2 ? (pos2 - pos1) == 1 : (pos1 - pos2) == 1;
  }

  Chain_barcode_swap& operator=(Chain_barcode_swap other) {
    Chain_pairing<Master_matrix>::operator=(other);
    pivotToPosition_.swap(other.pivotToPosition_);
  }
  friend void swap(Chain_barcode_swap& swap1, Chain_barcode_swap& swap2) {
    swap(static_cast<Chain_pairing<Master_matrix>&>(swap1), static_cast<Chain_pairing<Master_matrix>&>(swap2));
    swap1.pivotToPosition_.swap(swap2.pivotToPosition_);
  }

  Pos_index death(ID_index pivot) const {
    Pos_index simplexIndex = _get_pivot_position(pivot);

    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return CP::indexToBar_.at(simplexIndex)->death;
    } else {
      return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).death;
    }
  }

  Pos_index birth(ID_index pivot) const {
    Pos_index simplexIndex = _get_pivot_position(pivot);

    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return CP::indexToBar_.at(simplexIndex)->birth;
    } else {
      return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).birth;
    }
  }

 private:
  Pos_index _get_pivot_position(ID_index pivot) const {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      return pivotToPosition_.at(
          pivot);  // quite often called, make public and pass position instead of pivot to avoid find() every time?
    } else {
      return pivotToPosition_[pivot];
    }
  }

  Pos_index& _death(Pos_index simplexIndex) {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return CP::indexToBar_.at(simplexIndex)->death;
    } else {
      return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).death;
    }
  }

  Pos_index& _birth(Pos_index simplexIndex) {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return CP::indexToBar_.at(simplexIndex)->birth;
    } else {
      return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).birth;
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
class Chain_vine_swap : public std::conditional<Master_matrix::Option_list::has_column_pairings,
                                                Chain_barcode_swap<Master_matrix>, 
                                                Dummy_chain_vine_pairing
                                               >::type 
{
 public:
  using Index = typename Master_matrix::Index;                        /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;                  /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index;                /**< @ref PosIdx index type. */
  using Column_container = typename Master_matrix::Column_container;  /**< Column container type. */
  using Column = typename Master_matrix::Column;                      /**< Column type. */
  typedef bool (*EventCompFuncPointer)(Pos_index, Pos_index);         /**< Pointer type for birth/death comparators. */

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
  Chain_vine_swap(std::function<bool(Pos_index,Pos_index)> birthComparator,
                  std::function<bool(Pos_index,Pos_index)> deathComparator = _no_G_death_comparator);
  /**
   * @brief Copy constructor.
   * 
   * @param matrixToCopy Matrix to copy.
   */
  Chain_vine_swap(const Chain_vine_swap& matrixToCopy);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Chain_vine_swap(Chain_vine_swap&& other) noexcept;

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
   * @brief Assign operator.
   */
  Chain_vine_swap& operator=(Chain_vine_swap other);
  /**
   * @brief Swap operator.
   */
  friend void swap(Chain_vine_swap& swap1, Chain_vine_swap& swap2) {
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      swap(static_cast<Chain_barcode_swap<Master_matrix>&>(swap1),
           static_cast<Chain_barcode_swap<Master_matrix>&>(swap2));
    }
    std::swap(swap1.birthComp_, swap2.birthComp_);
    std::swap(swap1.deathComp_, swap2.deathComp_);
  }

 protected:
  using CP = typename std::conditional<Master_matrix::Option_list::has_column_pairings,
                                       Chain_barcode_swap<Master_matrix>, 
                                       Dummy_chain_vine_pairing
                                      >::type;

 private:
  using Master_chain_matrix = typename Master_matrix::Master_chain_matrix;

  std::function<bool(Pos_index,Pos_index)> birthComp_;  /**< for F x F & H x H. */
  std::function<bool(Pos_index,Pos_index)> deathComp_;  /**< for G x G. */

  bool _is_negative_in_pair(Index columnIndex);

  Index _positive_vine_swap(Index columnIndex1, Index columnIndex2);
  Index _positive_negative_vine_swap(Index columnIndex1, Index columnIndex2);
  Index _negative_positive_vine_swap(Index columnIndex1, Index columnIndex2);
  Index _negative_vine_swap(Index columnIndex1, Index columnIndex2);

  constexpr Master_chain_matrix* _matrix() { return static_cast<Master_chain_matrix*>(this); }
  constexpr const Master_chain_matrix* _matrix() const { return static_cast<const Master_chain_matrix*>(this); }
};

template <class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap() : CP(), birthComp_(), deathComp_() 
{
  static_assert(Master_matrix::Option_list::has_column_pairings,
                "If barcode is not stored, at least a birth comparator has to be specified.");
}

template <class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(std::function<bool(Pos_index,Pos_index)> birthComparator,
                                                       std::function<bool(Pos_index,Pos_index)> deathComparator)
    : CP(), birthComp_(std::move(birthComparator)), deathComp_(std::move(deathComparator)) 
{}

template <class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(const Chain_vine_swap& matrixToCopy)
    : CP(static_cast<const CP&>(matrixToCopy)),
      birthComp_(matrixToCopy.birthComp_),
      deathComp_(matrixToCopy.deathComp_) 
{}

template <class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(Chain_vine_swap<Master_matrix>&& other) noexcept
    : CP(std::move(static_cast<CP&>(other))),
      birthComp_(std::move(other.birthComp_)),
      deathComp_(std::move(other.deathComp_)) 
{}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(
    Index columnIndex1, Index columnIndex2) 
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
    GUDHI_CHECK(CP::are_adjacent(_matrix()->get_pivot(columnIndex1), _matrix()->get_pivot(columnIndex2)),
                std::invalid_argument(
                    "Chain_vine_swap::vine_swap - Columns to be swapped need to be adjacent in the 'real' matrix."));
  }

  const bool col1IsNeg = _is_negative_in_pair(columnIndex1);
  const bool col2IsNeg = _is_negative_in_pair(columnIndex2);

  if (col1IsNeg && col2IsNeg) {
    if (_matrix()->is_zero_entry(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
        ID_index pivot2 = _matrix()->get_pivot(columnIndex2);

        CP::negative_transpose(pivot1, pivot2);
        CP::swap_positions(pivot1, pivot2);
      }
      return columnIndex1;
    }
    return _negative_vine_swap(columnIndex1, columnIndex2);
  }

  if (col1IsNeg) {
    if (_matrix()->is_zero_entry(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
        ID_index pivot2 = _matrix()->get_pivot(columnIndex2);

        CP::negative_positive_transpose(pivot1, pivot2);
        CP::swap_positions(pivot1, pivot2);
      }
      return columnIndex1;
    }
    return _negative_positive_vine_swap(columnIndex1, columnIndex2);
  }

  if (col2IsNeg) {
    if (_matrix()->is_zero_entry(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
        ID_index pivot2 = _matrix()->get_pivot(columnIndex2);

        CP::positive_negative_transpose(pivot1, pivot2);
        CP::swap_positions(pivot1, pivot2);
      }
      return columnIndex1;
    }
    return _positive_negative_vine_swap(columnIndex1, columnIndex2);
  }

  if (_matrix()->is_zero_entry(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
      ID_index pivot2 = _matrix()->get_pivot(columnIndex2);

      CP::positive_transpose(pivot1, pivot2);
      CP::swap_positions(pivot1, pivot2);
    }
    return columnIndex1;
  }
  return _positive_vine_swap(columnIndex1, columnIndex2);
}

template <class Master_matrix>
inline Chain_vine_swap<Master_matrix>& Chain_vine_swap<Master_matrix>::operator=(Chain_vine_swap<Master_matrix> other) 
{
  CP::operator=(other);
  std::swap(birthComp_, other.birthComp_);
  std::swap(deathComp_, other.deathComp_);
  return *this;
}

template <class Master_matrix>
inline bool Chain_vine_swap<Master_matrix>::_is_negative_in_pair(Index columnIndex) 
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    return CP::is_negative_in_pair(_matrix()->get_pivot(columnIndex));
  } else {
    auto& col = _matrix()->get_column(columnIndex);
    if (!col.is_paired()) return false;
    return col.get_pivot() > _matrix()->get_pivot(col.get_paired_chain_index());
  }
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::_positive_vine_swap(
    Index columnIndex1, Index columnIndex2) 
{
  auto& col1 = _matrix()->get_column(columnIndex1);
  auto& col2 = _matrix()->get_column(columnIndex2);

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    CP::swap_positions(col1.get_pivot(), col2.get_pivot());
  }
  // TODO: factorize the cases. But for debug it is much more easier to understand what is happening splitted like this
  if (!col1.is_paired()) {  // F x *
    bool hasSmallerBirth;
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      // this order because position were swapped with CP::swap_positions
      hasSmallerBirth = (CP::birth(col2.get_pivot()) < CP::birth(col1.get_pivot()));
    } else {
      hasSmallerBirth = birthComp_(columnIndex1, columnIndex2);
    }

    if (!col2.is_paired() && hasSmallerBirth) {
      _matrix()->add_to(columnIndex1, columnIndex2);
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        CP::positive_transpose(col1.get_pivot(), col2.get_pivot());
      }
      return columnIndex1;
    }
    _matrix()->add_to(columnIndex2, columnIndex1);

    return columnIndex2;
  }

  if (!col2.is_paired()) {  // G x F
    static_cast<Master_chain_matrix*>(this)->add_to(columnIndex1, columnIndex2);
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      CP::positive_transpose(col1.get_pivot(), col2.get_pivot());
    }
    return columnIndex1;
  }

  bool hasSmallerDeath;
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    // this order because position were swapped with CP::swap_positions
    hasSmallerDeath = (CP::death(col2.get_pivot()) < CP::death(col1.get_pivot()));
  } else {
    hasSmallerDeath = deathComp_(columnIndex1, columnIndex2);
  }

  // G x G
  if (hasSmallerDeath)
  {
    _matrix()->add_to(col1.get_paired_chain_index(), col2.get_paired_chain_index());
    _matrix()->add_to(columnIndex1, columnIndex2);
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      CP::positive_transpose(col1.get_pivot(), col2.get_pivot());
    }
    return columnIndex1;
  }

  _matrix()->add_to(col2.get_paired_chain_index(), col1.get_paired_chain_index());
  _matrix()->add_to(columnIndex2, columnIndex1);

  return columnIndex2;
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::_positive_negative_vine_swap(
    Index columnIndex1, Index columnIndex2) 
{
  _matrix()->add_to(columnIndex1, columnIndex2);

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    ID_index pivot1 = _matrix()->get_pivot(columnIndex1);
    ID_index pivot2 = _matrix()->get_pivot(columnIndex2);

    CP::positive_negative_transpose(pivot1, pivot2);
    CP::swap_positions(pivot1, pivot2);
  }

  return columnIndex1;
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::_negative_positive_vine_swap(
    Index columnIndex1, Index columnIndex2) 
{
  _matrix()->add_to(columnIndex2, columnIndex1);

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    CP::swap_positions(_matrix()->get_pivot(columnIndex1), _matrix()->get_pivot(columnIndex2));
  }

  return columnIndex2;
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::Index Chain_vine_swap<Master_matrix>::_negative_vine_swap(
    Index columnIndex1, Index columnIndex2) 
{
  auto& col1 = _matrix()->get_column(columnIndex1);
  auto& col2 = _matrix()->get_column(columnIndex2);

  Index pairedIndex1 = col1.get_paired_chain_index();
  Index pairedIndex2 = col2.get_paired_chain_index();

  bool hasSmallerBirth;
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    hasSmallerBirth = (CP::birth(col1.get_pivot()) < CP::birth(col2.get_pivot()));
  } else {
    hasSmallerBirth = birthComp_(pairedIndex1, pairedIndex2);
  }

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    CP::swap_positions(col1.get_pivot(), col2.get_pivot());
  }

  if (hasSmallerBirth)
  {
    _matrix()->add_to(pairedIndex1, pairedIndex2);
    _matrix()->add_to(columnIndex1, columnIndex2);

    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      CP::negative_transpose(col1.get_pivot(), col2.get_pivot());
    }

    return columnIndex1;
  }

  _matrix()->add_to(pairedIndex2, pairedIndex1);
  _matrix()->add_to(columnIndex2, columnIndex1);

  return columnIndex2;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_CHAIN_VINE_SWAP_H
