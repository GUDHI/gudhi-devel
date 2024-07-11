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
 * @file ru_vine_swap.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::RU_vine_swap class, as well as the
 * @ref Gudhi::persistence_matrix::Dummy_ru_vine_swap and
 * @ref Gudhi::persistence_matrix::Dummy_ru_vine_pairing structures.
 */

#ifndef PM_RU_VINE_SWAP_H
#define PM_RU_VINE_SWAP_H

#include <utility>      //std::move
#include <type_traits>  //std::conditional
#include <cassert>
#include <vector>
#include <stdexcept>    //std::invalid_argument

#include "ru_pairing.h"

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inheritated instead of @ref RU_vine_swap, when vine swappes are not enabled.
 */
struct Dummy_ru_vine_swap {
  friend void swap([[maybe_unused]] Dummy_ru_vine_swap& d1, [[maybe_unused]] Dummy_ru_vine_swap& d2) {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inheritated instead of @ref RU_pairing, when the barcode is not stored.
 */
struct Dummy_ru_vine_pairing {
  friend void swap([[maybe_unused]] Dummy_ru_vine_pairing& d1, [[maybe_unused]] Dummy_ru_vine_pairing& d2) {}
};

/**
 * @class RU_vine_swap ru_vine_swap.h gudhi/Persistence_matrix/ru_vine_swap.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the vine swaps for @ref RU_matrix.
 * 
 * @tparam Master_matrix An instanciation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class RU_vine_swap : public std::conditional<Master_matrix::Option_list::has_column_pairings, 
                                             RU_pairing<Master_matrix>,
                                             Dummy_ru_vine_pairing
                                            >::type 
{
 public:
  using index = typename Master_matrix::index;          /**< @ref MatIdx index type. */
  using id_index = typename Master_matrix::id_index;    /**< @ref IDIdx index type. */
  using pos_index = typename Master_matrix::pos_index;  /**< @ref PosIdx index type. */

  /**
   * @brief Default constructor.
   */
  RU_vine_swap();
  /**
   * @brief Copy constructor.
   * 
   * @param matrixToCopy Matrix to copy.
   */
  RU_vine_swap(const RU_vine_swap& matrixToCopy);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  RU_vine_swap(RU_vine_swap&& other) noexcept;

  /**
   * @brief Does the same than @ref vine_swap, but assumes that the swap is non trivial and
   * therefore skips a part of the case study.
   * 
   * @param index @ref PosIdx index of the first face to swap. The second one has to be at `position + 1`.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap_with_z_eq_1_case(pos_index index);
  /**
   * @brief Does a vine swap between two faces which are consecutives in the filtration.
   * Roughly, if \f$ F \f$ is the current filtration represented by the matrix, the method modifies the matrix
   * such that the new state corresponds to a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but
   * with the two faces at position `position` and `position + 1` swapped. Of course, the two faces should
   * not have a face/coface relation which each other ; \f$ F' \f$ has to be a valid filtration.
   * See @cite vineyards for more information about vine and vineyards.
   * 
   * @param index @ref PosIdx index of the first face to swap. The second one has to be at `position + 1`.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap(pos_index index);

  /**
   * @brief Assign operator.
   */
  RU_vine_swap& operator=(RU_vine_swap other);
  /**
   * @brief Swap operator.
   */
  friend void swap(RU_vine_swap& swap1, RU_vine_swap& swap2) {
    if constexpr (Master_matrix::Option_list::has_column_pairings) {
      swap(static_cast<RU_pairing<Master_matrix>&>(swap1), static_cast<RU_pairing<Master_matrix>&>(swap2));
    }
    swap1.positionToRowIdx_.swap(swap2.positionToRowIdx_);
  }

 protected:
  // only usefull when simplex id does not corresponds to position, so feels kinda useless most of the time...
  // TODO: as it takes up some non trivial memory, see if this should not be optional
  // or only remember the positions with a difference. but then a map is needed, ie find instead of [].
  std::vector<id_index> positionToRowIdx_;  /**< Map from @ref PosIdx index to row index. */

 private:
  using RUP = typename std::conditional<Master_matrix::Option_list::has_column_pairings, 
                                        RU_pairing<Master_matrix>,
                                        Dummy_ru_vine_pairing
                                       >::type;
  using ru_matrix = typename Master_matrix::RU_matrix_type;

  bool _is_paired(index columnIndex);

  void _swap_at_index(index columnIndex);
  void _add_to(index sourceIndex, index targetIndex);
  void _positive_transpose(index columnIndex);
  void _negative_transpose(index columnIndex);
  void _positive_negative_transpose(index columnIndex);
  void _negative_positive_transpose(index columnIndex);
  bool _positive_vine_swap(index columnIndex);
  bool _negative_vine_swap(index columnIndex);
  bool _positive_negative_vine_swap(index columnIndex);
  bool _negative_positive_vine_swap(index columnIndex);

  pos_index& _death(pos_index simplexIndex);
  pos_index& _birth(pos_index simplexIndex);
  pos_index _get_death(index simplexIndex);
  pos_index _get_birth(index simplexIndex);

  constexpr ru_matrix* _matrix() { return static_cast<ru_matrix*>(this); }
  constexpr const ru_matrix* _matrix() const { return static_cast<const ru_matrix*>(this); }
};

template <class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap() : RUP() 
{}

template <class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(const RU_vine_swap& matrixToCopy)
    : RUP(static_cast<const RUP&>(matrixToCopy)), positionToRowIdx_(matrixToCopy.positionToRowIdx_) 
{}

template <class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(RU_vine_swap<Master_matrix>&& other) noexcept
    : RUP(std::move(static_cast<RUP&>(other))), positionToRowIdx_(std::move(other.positionToRowIdx_)) 
{}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(pos_index index) 
{
  GUDHI_CHECK(index < _matrix()->reducedMatrixR_.get_number_of_columns() - 1,
              std::invalid_argument("RU_vine_swap::vine_swap_with_z_eq_1_case - Index to swap out of bound."));

  bool iIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index);
  bool iiIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index + 1);

  if (iIsPositive && iiIsPositive) {
    _matrix()->mirrorMatrixU_.zero_cell(index, positionToRowIdx_[index + 1]);
    return _positive_vine_swap(index);
  } else if (!iIsPositive && !iiIsPositive) {
    return _negative_vine_swap(index);
  } else if (iIsPositive && !iiIsPositive) {
    return _positive_negative_vine_swap(index);
  } else {
    return _negative_positive_vine_swap(index);
  }
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::vine_swap(pos_index index) 
{
  GUDHI_CHECK(index < _matrix()->reducedMatrixR_.get_number_of_columns() - 1,
              std::invalid_argument("RU_vine_swap::vine_swap - Index to swap out of bound."));

  bool iIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index);
  bool iiIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index + 1);

  if (iIsPositive && iiIsPositive) {
    if (_matrix()->reducedMatrixR_.get_column_dimension(index) !=
        _matrix()->reducedMatrixR_.get_column_dimension(index + 1)) {
      _positive_transpose(index);
      _swap_at_index(index);
      return true;
    }
    if (!_matrix()->mirrorMatrixU_.is_zero_cell(index, positionToRowIdx_[index + 1])) {
      _matrix()->mirrorMatrixU_.zero_cell(index, positionToRowIdx_[index + 1]);
    }
    return _positive_vine_swap(index);
  } else if (!iIsPositive && !iiIsPositive) {
    if (_matrix()->reducedMatrixR_.get_column_dimension(index) !=
            _matrix()->reducedMatrixR_.get_column_dimension(index + 1) ||
        _matrix()->mirrorMatrixU_.is_zero_cell(index, positionToRowIdx_[index + 1])) {
      _negative_transpose(index);
      _swap_at_index(index);
      return true;
    }
    return _negative_vine_swap(index);
  } else if (iIsPositive && !iiIsPositive) {
    if (_matrix()->reducedMatrixR_.get_column_dimension(index) !=
            _matrix()->reducedMatrixR_.get_column_dimension(index + 1) ||
        _matrix()->mirrorMatrixU_.is_zero_cell(index, positionToRowIdx_[index + 1])) {
      _positive_negative_transpose(index);
      _swap_at_index(index);
      return true;
    }
    return _positive_negative_vine_swap(index);
  } else {
    if (_matrix()->reducedMatrixR_.get_column_dimension(index) !=
            _matrix()->reducedMatrixR_.get_column_dimension(index + 1) ||
        _matrix()->mirrorMatrixU_.is_zero_cell(index, positionToRowIdx_[index + 1])) {
      _negative_positive_transpose(index);
      _swap_at_index(index);
      return true;
    }
    return _negative_positive_vine_swap(index);
  }
}

template <class Master_matrix>
inline RU_vine_swap<Master_matrix>& RU_vine_swap<Master_matrix>::operator=(RU_vine_swap<Master_matrix> other) 
{
  RUP::operator=(other);
  positionToRowIdx_.swap(other.positionToRowIdx_);
  return *this;
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_is_paired(index columnIndex) 
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    return _get_death(columnIndex) != static_cast<pos_index>(-1);
  } else {
    if (!_matrix()->reducedMatrixR_.is_zero_column(columnIndex)) return true;

    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      if (_matrix()->pivotToColumnIndex_.find(columnIndex) == _matrix()->pivotToColumnIndex_.end()) return false;
    } else {
      if (_matrix()->pivotToColumnIndex_.operator[](columnIndex) == static_cast<index>(-1)) return false;
    }

    return true;
  }
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_swap_at_index(index columnIndex) 
{
  _matrix()->reducedMatrixR_.swap_columns(columnIndex, columnIndex + 1);
  _matrix()->reducedMatrixR_.swap_rows(positionToRowIdx_[columnIndex], positionToRowIdx_[columnIndex + 1]);
  _matrix()->mirrorMatrixU_.swap_columns(columnIndex, columnIndex + 1);
  _matrix()->mirrorMatrixU_.swap_rows(positionToRowIdx_[columnIndex], positionToRowIdx_[columnIndex + 1]);
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_add_to(index sourceIndex, index targetIndex) 
{
  _matrix()->reducedMatrixR_.add_to(sourceIndex, targetIndex);
  _matrix()->mirrorMatrixU_.add_to(targetIndex, sourceIndex);
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_transpose(index columnIndex) 
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if (_is_paired(columnIndex) && _is_paired(columnIndex + 1)) {
      std::swap(_matrix()->pivotToColumnIndex_.at(columnIndex), _matrix()->pivotToColumnIndex_.at(columnIndex + 1));
    } else if (_is_paired(columnIndex)) {
      _matrix()->pivotToColumnIndex_.emplace(columnIndex + 1, _matrix()->pivotToColumnIndex_.at(columnIndex));
      _matrix()->pivotToColumnIndex_.erase(columnIndex);
    } else if (_is_paired(columnIndex + 1)) {
      _matrix()->pivotToColumnIndex_.emplace(columnIndex, _matrix()->pivotToColumnIndex_.at(columnIndex + 1));
      _matrix()->pivotToColumnIndex_.erase(columnIndex + 1);
    }
  } else {
    std::swap(_matrix()->pivotToColumnIndex_.operator[](columnIndex),
              _matrix()->pivotToColumnIndex_.operator[](columnIndex + 1));
  }

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _birth(columnIndex) = columnIndex + 1;
    _birth(columnIndex + 1) = columnIndex;
    std::swap(RUP::indexToBar_.at(columnIndex), RUP::indexToBar_.at(columnIndex + 1));
  }
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_transpose(index columnIndex) 
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _death(columnIndex) = columnIndex + 1;
    _death(columnIndex + 1) = columnIndex;
    std::swap(RUP::indexToBar_.at(columnIndex), RUP::indexToBar_.at(columnIndex + 1));
  }
  std::swap(_matrix()->pivotToColumnIndex_.at(_get_birth(columnIndex)),
            _matrix()->pivotToColumnIndex_.at(_get_birth(columnIndex + 1)));
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_negative_transpose(index columnIndex) 
{
  _matrix()->pivotToColumnIndex_.at(_get_birth(columnIndex + 1)) = columnIndex;
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if (_is_paired(columnIndex)) {
      _matrix()->pivotToColumnIndex_.emplace(columnIndex + 1, _matrix()->pivotToColumnIndex_.at(columnIndex));
      _matrix()->pivotToColumnIndex_.erase(columnIndex);
    }
  } else {
    _matrix()->pivotToColumnIndex_.operator[](columnIndex + 1) = _matrix()->pivotToColumnIndex_.operator[](columnIndex);
    _matrix()->pivotToColumnIndex_.operator[](columnIndex) = -1;
  }

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _birth(columnIndex) = columnIndex + 1;
    _death(columnIndex + 1) = columnIndex;
    std::swap(RUP::indexToBar_.at(columnIndex), RUP::indexToBar_.at(columnIndex + 1));
  }
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_positive_transpose(index columnIndex) 
{
  _matrix()->pivotToColumnIndex_.at(_get_birth(columnIndex)) = columnIndex + 1;
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    if (_is_paired(columnIndex + 1)) {
      _matrix()->pivotToColumnIndex_.emplace(columnIndex, _matrix()->pivotToColumnIndex_.at(columnIndex + 1));
      _matrix()->pivotToColumnIndex_.erase(columnIndex + 1);
    }
  } else {
    _matrix()->pivotToColumnIndex_.operator[](columnIndex) = _matrix()->pivotToColumnIndex_.operator[](columnIndex + 1);
    _matrix()->pivotToColumnIndex_.operator[](columnIndex + 1) = -1;
  }

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _death(columnIndex) = columnIndex + 1;
    _birth(columnIndex + 1) = columnIndex;
    std::swap(RUP::indexToBar_.at(columnIndex), RUP::indexToBar_.at(columnIndex + 1));
  }
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_positive_vine_swap(index columnIndex) 
{
  const pos_index iDeath = _get_death(columnIndex);
  const pos_index iiDeath = _get_death(columnIndex + 1);

  if (iDeath != static_cast<pos_index>(-1) && iiDeath != static_cast<pos_index>(-1) &&
      !(_matrix()->reducedMatrixR_.is_zero_cell(iiDeath, positionToRowIdx_[columnIndex]))) {
    if (iDeath < iiDeath) {
      _swap_at_index(columnIndex);
      _add_to(iDeath, iiDeath);
      _positive_transpose(columnIndex);
      return true;
    }

    _swap_at_index(columnIndex);
    _add_to(iiDeath, iDeath);
    return false;
  }

  _swap_at_index(columnIndex);

  if (iDeath != static_cast<pos_index>(-1) || iiDeath == static_cast<pos_index>(-1) ||
      _matrix()->reducedMatrixR_.is_zero_cell(iiDeath, positionToRowIdx_[columnIndex + 1])) {
    _positive_transpose(columnIndex);
    return true;
  }
  return false;
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_negative_vine_swap(index columnIndex) 
{
  const pos_index iBirth = _get_birth(columnIndex);
  const pos_index iiBirth = _get_birth(columnIndex + 1);

  _add_to(columnIndex, columnIndex + 1);
  _swap_at_index(columnIndex);

  if (iBirth < iiBirth) {
    _negative_transpose(columnIndex);
    return true;
  }

  _add_to(columnIndex, columnIndex + 1);

  return false;
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_positive_negative_vine_swap(index columnIndex) 
{
  _matrix()->mirrorMatrixU_.zero_cell(columnIndex, positionToRowIdx_[columnIndex + 1]);

  _swap_at_index(columnIndex);
  _positive_negative_transpose(columnIndex);

  return true;
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_negative_positive_vine_swap(index columnIndex) 
{
  _add_to(columnIndex, columnIndex + 1);  // useless for R?
  _swap_at_index(columnIndex);      // if additions not made for R, do not swap R columns, just rows
  _add_to(columnIndex, columnIndex + 1);  // useless for R?

  return false;
}

template <class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::pos_index& RU_vine_swap<Master_matrix>::_death(pos_index simplexIndex) 
{
  static_assert(Master_matrix::Option_list::has_column_pairings, "Pairing necessary to modify death value.");

  if constexpr (Master_matrix::Option_list::has_removable_columns) {
    return RUP::indexToBar_.at(simplexIndex)->death;
  } else {
    return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).death;
  }
}

template <class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::pos_index& RU_vine_swap<Master_matrix>::_birth(pos_index simplexIndex) 
{
  static_assert(Master_matrix::Option_list::has_column_pairings, "Pairing necessary to modify birth value.");

  if constexpr (Master_matrix::Option_list::has_removable_columns) {
    return RUP::indexToBar_.at(simplexIndex)->birth;
  } else {
    return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).birth;
  }
}

template <class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::pos_index RU_vine_swap<Master_matrix>::_get_death(index simplexIndex) 
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return RUP::indexToBar_.at(simplexIndex)->death;
    } else {
      return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).death;
    }
  } else {
    if (!_matrix()->reducedMatrixR_.is_zero_column(simplexIndex))
      return _matrix()->reducedMatrixR_.get_column(simplexIndex).get_pivot();

    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      auto it = _matrix()->pivotToColumnIndex_.find(simplexIndex);
      if (it == _matrix()->pivotToColumnIndex_.end()) return -1;
      return it->second;
    } else {
      return _matrix()->pivotToColumnIndex_.operator[](simplexIndex);
    }
  }
}

template <class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::pos_index RU_vine_swap<Master_matrix>::_get_birth(
    index negativeSimplexIndex) 
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return RUP::indexToBar_.at(negativeSimplexIndex)->birth;
    } else {
      return RUP::barcode_.at(RUP::indexToBar_.at(negativeSimplexIndex)).birth;
    }
  } else {
    return _matrix()->reducedMatrixR_.get_pivot(negativeSimplexIndex);
  }
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_RU_VINE_SWAP_H
