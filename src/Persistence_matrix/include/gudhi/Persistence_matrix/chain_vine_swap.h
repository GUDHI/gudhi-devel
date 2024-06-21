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
 * @brief Contains the @ref Chain_barcode_swap and @ref Chain_vine_swap classes, as well as the
 * @ref Dummy_chain_vine_swap and @ref Dummy_chain_vine_pairing structures.
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
 * Inheritated instead of @ref Chain_vine_swap, when vine swappes are not enabled.
 */
struct Dummy_chain_vine_swap {
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
 * Inheritated instead of @ref Chain_barcode_swap, when the barcode is not stored.
 */
struct Dummy_chain_vine_pairing {
  friend void swap([[maybe_unused]] Dummy_chain_vine_pairing& d1, [[maybe_unused]] Dummy_chain_vine_pairing& d2) {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Class managing the barcode for @ref Chain_vine_swap.
 * 
 * @tparam Master_matrix An instanciation of @ref Matrix from which all types and options are deduced.
 */
template <typename Master_matrix>
class Chain_barcode_swap : public Chain_pairing<Master_matrix> 
{
 public:
  using id_index = typename Master_matrix::id_index;    /**< @ref IDIdx index type. */
  using pos_index = typename Master_matrix::pos_index;  /**< @ref PosIdx index type. */
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
  using dictionnary_type = typename Master_matrix::template dictionnary_type<pos_index>;

  dictionnary_type pivotToPosition_;  // necessary to keep track of the barcode changes

  void swap_positions(id_index pivot1, id_index pivot2) {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      std::swap(pivotToPosition_.at(pivot1), pivotToPosition_.at(pivot2));
    } else {
      std::swap(pivotToPosition_[pivot1], pivotToPosition_[pivot2]);
    }
  }

  bool is_negative_in_pair(id_index pivot) const {
    pos_index pos = _get_pivot_position(pivot);
    return death(pivot) == pos;
  }

  void positive_transpose(id_index pivot1, id_index pivot2) {
    pos_index pos1 = _get_pivot_position(pivot1);
    pos_index pos2 = _get_pivot_position(pivot2);

    _birth(pos1) = pos2;
    _birth(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  void negative_transpose(id_index pivot1, id_index pivot2) {
    pos_index pos1 = _get_pivot_position(pivot1);
    pos_index pos2 = _get_pivot_position(pivot2);

    _death(pos1) = pos2;
    _death(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  void positive_negative_transpose(id_index pivot1, id_index pivot2) {
    pos_index pos1 = _get_pivot_position(pivot1);
    pos_index pos2 = _get_pivot_position(pivot2);

    _birth(pos1) = pos2;
    _death(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  void negative_positive_transpose(id_index pivot1, id_index pivot2) {
    pos_index pos1 = _get_pivot_position(pivot1);
    pos_index pos2 = _get_pivot_position(pivot2);

    _death(pos1) = pos2;
    _birth(pos2) = pos1;
    std::swap(CP::indexToBar_.at(pos1), CP::indexToBar_.at(pos2));
  }

  bool are_adjacent(id_index pivot1, id_index pivot2) const {
    pos_index pos1 = _get_pivot_position(pivot1);
    pos_index pos2 = _get_pivot_position(pivot2);
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

  pos_index death(id_index pivot) const {
    pos_index simplexIndex = _get_pivot_position(pivot);

    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return CP::indexToBar_.at(simplexIndex)->death;
    } else {
      return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).death;
    }
  }

  pos_index birth(id_index pivot) const {
    pos_index simplexIndex = _get_pivot_position(pivot);

    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return CP::indexToBar_.at(simplexIndex)->birth;
    } else {
      return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).birth;
    }
  }

 private:
  pos_index _get_pivot_position(id_index pivot) const {
    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      return pivotToPosition_.at(
          pivot);  // quite often called, make public and pass position instead of pivot to avoid find() everytime?
    } else {
      return pivotToPosition_[pivot];
    }
  }

  pos_index& _death(pos_index simplexIndex) {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return CP::indexToBar_.at(simplexIndex)->death;
    } else {
      return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).death;
    }
  }

  pos_index& _birth(pos_index simplexIndex) {
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
 * @tparam Master_matrix An instanciation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Chain_vine_swap : public std::conditional<Master_matrix::Option_list::has_column_pairings,
                                                Chain_barcode_swap<Master_matrix>, 
                                                Dummy_chain_vine_pairing
                                               >::type 
{
 public:
  using index = typename Master_matrix::index;                        /**< @ref MatIdx index type. */
  using id_index = typename Master_matrix::id_index;                  /**< @ref IDIdx index type. */
  using pos_index = typename Master_matrix::pos_index;                /**< @ref PosIdx index type. */
  using matrix_type = typename Master_matrix::column_container_type;  /**< Column container type. */
  using Column_type = typename Master_matrix::Column_type;            /**< Column type. */
  typedef bool (*EventCompFuncPointer)(pos_index, pos_index);         /**< Pointer type for birth/death comparators. */

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
  Chain_vine_swap(std::function<bool(pos_index,pos_index)> birthComparator,
                  std::function<bool(pos_index,pos_index)> deathComparator = _no_G_death_comparator);
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
   * @param columnIndex1 @ref MatIdx index of the first face.
   * @param columnIndex2 @ref MatIdx index of the second face.
   * @return Let \f$ pos1 \f$ be the @ref PosIdx index of @p columnIndex1 and \f$ pos2 \f$ be the @ref PosIdx index of
   * @p columnIndex2. The method returns the @ref MatIdx of the column which has now, after the swap, the @ref PosIdx
   * \f$ max(pos1, pos2) \f$.
   */
  index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);
  /**
   * @brief Does a vine swap between two faces which are consecutives in the filtration. Roughly, if \f$ F \f$ is
   * the current filtration represented by the matrix, the method modifies the matrix such that the new state
   * corresponds to a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but with the two given faces
   * at swapped positions. Of course, the two faces should not have a face/coface relation which each other ;
   * \f$ F' \f$ has to be a valid filtration.
   * See @cite vineyards for more information about vine and vineyards.
   * 
   * @param columnIndex1 @ref MatIdx index of the first face.
   * @param columnIndex2 @ref MatIdx index of the second face. It is assumed that the @ref PosIdx of both only differs
   * by one if the barcode is maintained.
   * @return Let \f$ pos1 \f$ be the @ref PosIdx index of @p columnIndex1 and \f$ pos2 \f$ be the @ref PosIdx index of
   * @p columnIndex2. The method returns the @ref MatIdx of the column which has now, after the swap, the @ref PosIdx
   * \f$ max(pos1, pos2) \f$.
   */
  index vine_swap(index columnIndex1, index columnIndex2);

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
  using chain_matrix = typename Master_matrix::Chain_matrix_type;

  std::function<bool(pos_index,pos_index)> birthComp_;  /**< for F x F & H x H. */
  std::function<bool(pos_index,pos_index)> deathComp_;  /**< for G x G. */

  bool _is_negative_in_pair(index columnIndex);

  index _positive_vine_swap(index columnIndex1, index columnIndex2);
  index _positive_negative_vine_swap(index columnIndex1, index columnIndex2);
  index _negative_positive_vine_swap(index columnIndex1, index columnIndex2);
  index _negative_vine_swap(index columnIndex1, index columnIndex2);

  constexpr chain_matrix* _matrix() { return static_cast<chain_matrix*>(this); }
  constexpr const chain_matrix* _matrix() const { return static_cast<const chain_matrix*>(this); }
};

template <class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap() : CP(), birthComp_(), deathComp_() 
{
  static_assert(Master_matrix::Option_list::has_column_pairings,
                "If barcode is not stored, at least a birth comparator has to be specified.");
}

template <class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(std::function<bool(pos_index,pos_index)> birthComparator,
                                                       std::function<bool(pos_index,pos_index)> deathComparator)
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
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(
    index columnIndex1, index columnIndex2) 
{
  const bool col1IsNeg = _is_negative_in_pair(columnIndex1);
  const bool col2IsNeg = _is_negative_in_pair(columnIndex2);

  if (col1IsNeg && col2IsNeg) return _negative_vine_swap(columnIndex1, columnIndex2);

  if (col1IsNeg) return _negative_positive_vine_swap(columnIndex1, columnIndex2);

  if (col2IsNeg) return _positive_negative_vine_swap(columnIndex1, columnIndex2);

  return _positive_vine_swap(columnIndex1, columnIndex2);
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::vine_swap(index columnIndex1,
                                                                                                index columnIndex2) 
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    GUDHI_CHECK(CP::are_adjacent(_matrix()->get_pivot(columnIndex1), _matrix()->get_pivot(columnIndex2)),
                std::invalid_argument(
                    "Chain_vine_swap::vine_swap - Columns to be swaped need to be adjacent in the 'real' matrix."));
  }

  const bool col1IsNeg = _is_negative_in_pair(columnIndex1);
  const bool col2IsNeg = _is_negative_in_pair(columnIndex2);

  if (col1IsNeg && col2IsNeg) {
    if (_matrix()->is_zero_cell(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        id_index pivot1 = _matrix()->get_pivot(columnIndex1);
        id_index pivot2 = _matrix()->get_pivot(columnIndex2);

        CP::negative_transpose(pivot1, pivot2);
        CP::swap_positions(pivot1, pivot2);
      }
      return columnIndex1;
    }
    return _negative_vine_swap(columnIndex1, columnIndex2);
  }

  if (col1IsNeg) {
    if (_matrix()->is_zero_cell(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        id_index pivot1 = _matrix()->get_pivot(columnIndex1);
        id_index pivot2 = _matrix()->get_pivot(columnIndex2);

        CP::negative_positive_transpose(pivot1, pivot2);
        CP::swap_positions(pivot1, pivot2);
      }
      return columnIndex1;
    }
    return _negative_positive_vine_swap(columnIndex1, columnIndex2);
  }

  if (col2IsNeg) {
    if (_matrix()->is_zero_cell(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
      if constexpr (Master_matrix::Option_list::has_column_pairings) {
        id_index pivot1 = _matrix()->get_pivot(columnIndex1);
        id_index pivot2 = _matrix()->get_pivot(columnIndex2);

        CP::positive_negative_transpose(pivot1, pivot2);
        CP::swap_positions(pivot1, pivot2);
      }
      return columnIndex1;
    }
    return _positive_negative_vine_swap(columnIndex1, columnIndex2);
  }

  if (_matrix()->is_zero_cell(columnIndex2, _matrix()->get_pivot(columnIndex1))) {
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      id_index pivot1 = _matrix()->get_pivot(columnIndex1);
      id_index pivot2 = _matrix()->get_pivot(columnIndex2);

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
inline bool Chain_vine_swap<Master_matrix>::_is_negative_in_pair(index columnIndex) 
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
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::_positive_vine_swap(
    index columnIndex1, index columnIndex2) 
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
    static_cast<chain_matrix*>(this)->add_to(columnIndex1, columnIndex2);
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
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::_positive_negative_vine_swap(
    index columnIndex1, index columnIndex2) 
{
  _matrix()->add_to(columnIndex1, columnIndex2);

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    id_index pivot1 = _matrix()->get_pivot(columnIndex1);
    id_index pivot2 = _matrix()->get_pivot(columnIndex2);

    CP::positive_negative_transpose(pivot1, pivot2);
    CP::swap_positions(pivot1, pivot2);
  }

  return columnIndex1;
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::_negative_positive_vine_swap(
    index columnIndex1, index columnIndex2) 
{
  _matrix()->add_to(columnIndex2, columnIndex1);

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    CP::swap_positions(_matrix()->get_pivot(columnIndex1), _matrix()->get_pivot(columnIndex2));
  }

  return columnIndex2;
}

template <class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::_negative_vine_swap(
    index columnIndex1, index columnIndex2) 
{
  auto& col1 = _matrix()->get_column(columnIndex1);
  auto& col2 = _matrix()->get_column(columnIndex2);

  index pairedIndex1 = col1.get_paired_chain_index();
  index pairedIndex2 = col2.get_paired_chain_index();

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
