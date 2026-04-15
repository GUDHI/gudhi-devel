/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
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

#include <cassert>
#include <utility>    //std::move
#include <stdexcept>  //std::invalid_argument

#include <gudhi/Persistence_matrix/ru_pairing.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref RU_vine_swap, when vine swaps are not enabled.
 */
struct Dummy_ru_vine_swap {
  using RUP = void;

  friend void swap([[maybe_unused]] Dummy_ru_vine_swap& d1, [[maybe_unused]] Dummy_ru_vine_swap& d2) noexcept {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref RU_pairing, when the barcode is not stored.
 */
struct Dummy_ru_vine_pairing {
  friend void swap([[maybe_unused]] Dummy_ru_vine_pairing& d1, [[maybe_unused]] Dummy_ru_vine_pairing& d2) noexcept {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Class managing the barcode for @ref RU_vine_swap.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <typename Master_matrix>
class RU_barcode_swap : public RU_pairing<Master_matrix>
{
 public:
  using Index = typename Master_matrix::Index;         /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;   /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index; /**< @ref PosIdx index type. */
  // RUP = RU Pairing
  using RUP = RU_pairing<Master_matrix>;

  /**
   * @brief Default constructor.
   */
  RU_barcode_swap() = default;
  /**
   * @brief Copy constructor.
   *
   * @param toCopy Matrix to copy.
   */
  RU_barcode_swap(const RU_barcode_swap& toCopy) : RUP(static_cast<const RUP&>(toCopy)) {};
  /**
   * @brief Move constructor.
   *
   * @param other Matrix to move.
   */
  RU_barcode_swap(RU_barcode_swap&& other) noexcept : RUP(std::move(static_cast<RUP&>(other))) {};

  ~RU_barcode_swap() = default;

  RU_barcode_swap& operator=(const RU_barcode_swap& other)
  {
    RUP::operator=(other);
    return *this;
  }

  RU_barcode_swap& operator=(RU_barcode_swap&& other) noexcept
  {
    RUP::operator=(std::move(other));
    return *this;
  }

  friend void swap(RU_barcode_swap& swap1, RU_barcode_swap& swap2) noexcept
  {
    swap(static_cast<RUP&>(swap1), static_cast<RUP&>(swap2));
  }

 protected:
  void _positive_transpose_barcode(Index columnIndex)
  {
    _birth(columnIndex) = columnIndex + 1;
    _birth(columnIndex + 1) = columnIndex;
    std::swap(RUP::indexToBar_.at(columnIndex), RUP::indexToBar_.at(columnIndex + 1));
  }

  void _negative_transpose_barcode(Index columnIndex)
  {
    _death(columnIndex) = columnIndex + 1;
    _death(columnIndex + 1) = columnIndex;
    std::swap(RUP::indexToBar_.at(columnIndex), RUP::indexToBar_.at(columnIndex + 1));
  }

  void _positive_negative_transpose_barcode(Index columnIndex)
  {
    _birth(columnIndex) = columnIndex + 1;
    _death(columnIndex + 1) = columnIndex;
    std::swap(RUP::indexToBar_.at(columnIndex), RUP::indexToBar_.at(columnIndex + 1));
  }

  void _negative_positive_transpose_barcode(Index columnIndex)
  {
    _death(columnIndex) = columnIndex + 1;
    _birth(columnIndex + 1) = columnIndex;
    std::swap(RUP::indexToBar_.at(columnIndex), RUP::indexToBar_.at(columnIndex + 1));
  }

  Pos_index _death_val(Pos_index index) const
  {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return RUP::indexToBar_.at(index)->death;
    } else {
      return RUP::barcode_.at(RUP::indexToBar_.at(index)).death;
    }
  }

  Pos_index _birth_val(Pos_index index) const
  {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return RUP::indexToBar_.at(index)->birth;
    } else {
      return RUP::barcode_.at(RUP::indexToBar_.at(index)).birth;
    }
  }

  void _reset() { RUP::_reset(); }

 private:
  Pos_index& _death(Pos_index index)
  {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return RUP::indexToBar_.at(index)->death;
    } else {
      return RUP::barcode_.at(RUP::indexToBar_.at(index)).death;
    }
  }

  Pos_index& _birth(Pos_index index)
  {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      return RUP::indexToBar_.at(index)->birth;
    } else {
      return RUP::barcode_.at(RUP::indexToBar_.at(index)).birth;
    }
  }
};

/**
 * @class RU_vine_swap ru_vine_swap.h gudhi/Persistence_matrix/ru_vine_swap.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the vine swaps for @ref RU_matrix.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class RU_vine_swap
{
 public:
  using Index = typename Master_matrix::Index;         /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;   /**< @ref IDIdx index type. */
  using Pos_index = typename Master_matrix::Pos_index; /**< @ref PosIdx index type. */

  /**
   * @brief Default constructor.
   */
  RU_vine_swap() = default;

  /**
   * @brief Does the same than @ref vine_swap, but assumes that the swap is non trivial and
   * therefore skips a part of the case study.
   *
   * @param index @ref PosIdx index of the first cell to swap. The second one has to be at `position + 1`.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap_with_z_eq_1_case(Pos_index index);
  /**
   * @brief Does a vine swap between two cells which are consecutive in the filtration.
   * Roughly, if \f$ F \f$ is the current filtration represented by the matrix, the method modifies the matrix
   * such that the new state corresponds to a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but
   * with the two cells at position `position` and `position + 1` swapped. Of course, the two cells should
   * not have a face/coface relation which each other ; \f$ F' \f$ has to be a valid filtration.
   * See @cite vineyards for more information about vine and vineyards.
   *
   * @param index @ref PosIdx index of the first cell to swap. The second one has to be at `position + 1`.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap(Pos_index index);

  /**
   * @brief Swap operator.
   */
  friend void swap([[maybe_unused]] RU_vine_swap& swap1, [[maybe_unused]] RU_vine_swap& swap2) noexcept {}

 private:
  using Master_RU_matrix = typename Master_matrix::Master_RU_matrix;

  bool _is_paired(Index columnIndex);
  void _swap_at_index(Index columnIndex);
  void _add_to(Index sourceIndex, Index targetIndex);
  void _positive_transpose(Index columnIndex);
  void _negative_transpose(Index columnIndex);
  void _positive_negative_transpose(Index columnIndex);
  void _negative_positive_transpose(Index columnIndex);
  bool _positive_vine_swap(Index columnIndex);
  bool _negative_vine_swap(Index columnIndex);
  bool _positive_negative_vine_swap(Index columnIndex);
  bool _negative_positive_vine_swap(Index columnIndex);
  Pos_index _get_death(Index simplexIndex);
  Pos_index _get_birth(Index simplexIndex);
  ID_index _get_row_id_from_position(Pos_index position) const;

  constexpr Master_RU_matrix* _matrix() { return static_cast<Master_RU_matrix*>(this); }

  constexpr const Master_RU_matrix* _matrix() const { return static_cast<const Master_RU_matrix*>(this); }
};

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(Pos_index index)
{
  GUDHI_CHECK(index < _matrix()->reducedMatrixR_.get_number_of_columns() - 1,
              std::invalid_argument("RU_vine_swap::vine_swap_with_z_eq_1_case - Index to swap out of bound."));

  bool iIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index);
  bool iiIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index + 1);

  if (iIsPositive && iiIsPositive) {
    _matrix()->mirrorMatrixU_.zero_entry(index, _get_row_id_from_position(index + 1));
    return _positive_vine_swap(index);
  }
  if (!iIsPositive && !iiIsPositive) {
    return _negative_vine_swap(index);
  }
  if (iIsPositive && !iiIsPositive) {
    return _positive_negative_vine_swap(index);
  }
  return _negative_positive_vine_swap(index);
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::vine_swap(Pos_index index)
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
    if (!_matrix()->mirrorMatrixU_.is_zero_entry(index, _get_row_id_from_position(index + 1))) {
      _matrix()->mirrorMatrixU_.zero_entry(index, _get_row_id_from_position(index + 1));
    }
    return _positive_vine_swap(index);
  }
  if (!iIsPositive && !iiIsPositive) {
    if (_matrix()->reducedMatrixR_.get_column_dimension(index) !=
            _matrix()->reducedMatrixR_.get_column_dimension(index + 1) ||
        _matrix()->mirrorMatrixU_.is_zero_entry(index, _get_row_id_from_position(index + 1))) {
      _negative_transpose(index);
      _swap_at_index(index);
      return true;
    }
    return _negative_vine_swap(index);
  }
  if (iIsPositive && !iiIsPositive) {
    if (_matrix()->reducedMatrixR_.get_column_dimension(index) !=
            _matrix()->reducedMatrixR_.get_column_dimension(index + 1) ||
        _matrix()->mirrorMatrixU_.is_zero_entry(index, _get_row_id_from_position(index + 1))) {
      _positive_negative_transpose(index);
      _swap_at_index(index);
      return true;
    }
    return _positive_negative_vine_swap(index);
  }
  if (_matrix()->reducedMatrixR_.get_column_dimension(index) !=
          _matrix()->reducedMatrixR_.get_column_dimension(index + 1) ||
      _matrix()->mirrorMatrixU_.is_zero_entry(index, _get_row_id_from_position(index + 1))) {
    _negative_positive_transpose(index);
    _swap_at_index(index);
    return true;
  }
  return _negative_positive_vine_swap(index);
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_is_paired(Index columnIndex)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    return _get_death(columnIndex) != Master_matrix::template get_null_value<Pos_index>();
  } else {
    if (!_matrix()->reducedMatrixR_.is_zero_column(columnIndex)) return true;

    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      return _matrix()->pivotToColumnIndex_.find(columnIndex) != _matrix()->pivotToColumnIndex_.end();
    } else {
      if (_matrix()->pivotToColumnIndex_.size() <= columnIndex) return false;
      return _matrix()->pivotToColumnIndex_[columnIndex] != Master_matrix::template get_null_value<Index>();
    }
  }
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_swap_at_index(Index columnIndex)
{
  _matrix()->reducedMatrixR_.swap_columns(columnIndex, columnIndex + 1);
  _matrix()->reducedMatrixR_.swap_rows(_get_row_id_from_position(columnIndex),
                                       _get_row_id_from_position(columnIndex + 1));
  _matrix()->mirrorMatrixU_.swap_columns(columnIndex, columnIndex + 1);
  _matrix()->mirrorMatrixU_.swap_rows(_get_row_id_from_position(columnIndex),
                                      _get_row_id_from_position(columnIndex + 1));
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_add_to(Index sourceIndex, Index targetIndex)
{
  _matrix()->reducedMatrixR_.add_to(sourceIndex, targetIndex);
  _matrix()->mirrorMatrixU_.add_to(targetIndex, sourceIndex);
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_transpose(Index columnIndex)
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
    if (_is_paired(columnIndex) || _is_paired(columnIndex + 1)) {
      if (columnIndex + 1 >= _matrix()->pivotToColumnIndex_.size()) {
        // pivotToColumnIndex_ has at least size columnIndex + 1
        _matrix()->pivotToColumnIndex_.push_back(Master_matrix::template get_null_value<Index>());
      }
      std::swap(_matrix()->pivotToColumnIndex_[columnIndex], _matrix()->pivotToColumnIndex_[columnIndex + 1]);
    }
  }

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _matrix()->_positive_transpose_barcode(columnIndex);
  }
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_transpose(Index columnIndex)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _matrix()->_negative_transpose_barcode(columnIndex);
  }
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    std::swap(_matrix()->pivotToColumnIndex_.at(_get_birth(columnIndex)),
              _matrix()->pivotToColumnIndex_.at(_get_birth(columnIndex + 1)));
  } else {
    std::swap(_matrix()->pivotToColumnIndex_[_get_birth(columnIndex)],
              _matrix()->pivotToColumnIndex_[_get_birth(columnIndex + 1)]);
  }
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_negative_transpose(Index columnIndex)
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    _matrix()->pivotToColumnIndex_.at(_get_birth(columnIndex + 1)) = columnIndex;
    if (_is_paired(columnIndex)) {
      _matrix()->pivotToColumnIndex_.emplace(columnIndex + 1, _matrix()->pivotToColumnIndex_.at(columnIndex));
      _matrix()->pivotToColumnIndex_.erase(columnIndex);
    }
  } else {
    _matrix()->pivotToColumnIndex_[_get_birth(columnIndex + 1)] = columnIndex;
    if (_is_paired(columnIndex)){
      if (columnIndex + 1 >= _matrix()->pivotToColumnIndex_.size()) {
        // pivotToColumnIndex_ has at least size columnIndex + 1
        _matrix()->pivotToColumnIndex_.push_back(Master_matrix::template get_null_value<Index>());
      }
      _matrix()->pivotToColumnIndex_[columnIndex + 1] = _matrix()->pivotToColumnIndex_[columnIndex];
      _matrix()->pivotToColumnIndex_[columnIndex] = Master_matrix::template get_null_value<Index>();
    }
  }

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _matrix()->_positive_negative_transpose_barcode(columnIndex);
  }
}

template <class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_positive_transpose(Index columnIndex)
{
  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    _matrix()->pivotToColumnIndex_.at(_get_birth(columnIndex)) = columnIndex + 1;
    if (_is_paired(columnIndex + 1)) {
      _matrix()->pivotToColumnIndex_.emplace(columnIndex, _matrix()->pivotToColumnIndex_.at(columnIndex + 1));
      _matrix()->pivotToColumnIndex_.erase(columnIndex + 1);
    }
  } else {
    _matrix()->pivotToColumnIndex_[_get_birth(columnIndex)] = columnIndex + 1;
    if (_is_paired(columnIndex + 1)){
      _matrix()->pivotToColumnIndex_[columnIndex] = _matrix()->pivotToColumnIndex_[columnIndex + 1];
      _matrix()->pivotToColumnIndex_[columnIndex + 1] = Master_matrix::template get_null_value<Index>();
    }
  }

  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    _matrix()->_negative_positive_transpose_barcode(columnIndex);
  }
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_positive_vine_swap(Index columnIndex)
{
  const Pos_index iDeath = _get_death(columnIndex);
  const Pos_index iiDeath = _get_death(columnIndex + 1);

  if (iDeath != Master_matrix::template get_null_value<Pos_index>() &&
      iiDeath != Master_matrix::template get_null_value<Pos_index>() &&
      !(_matrix()->reducedMatrixR_.is_zero_entry(iiDeath, _get_row_id_from_position(columnIndex)))) {
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

  if (iDeath != Master_matrix::template get_null_value<Pos_index>() ||
      iiDeath == Master_matrix::template get_null_value<Pos_index>() ||
      _matrix()->reducedMatrixR_.is_zero_entry(iiDeath, _get_row_id_from_position(columnIndex + 1))) {
    _positive_transpose(columnIndex);
    return true;
  }
  return false;
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_negative_vine_swap(Index columnIndex)
{
  const Pos_index iBirth = _get_birth(columnIndex);
  const Pos_index iiBirth = _get_birth(columnIndex + 1);

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
inline bool RU_vine_swap<Master_matrix>::_positive_negative_vine_swap(Index columnIndex)
{
  _matrix()->mirrorMatrixU_.zero_entry(columnIndex, _get_row_id_from_position(columnIndex + 1));

  _positive_negative_transpose(columnIndex);
  _swap_at_index(columnIndex);

  return true;
}

template <class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_negative_positive_vine_swap(Index columnIndex)
{
  _add_to(columnIndex, columnIndex + 1);  // useless for R?
  _swap_at_index(columnIndex);            // if additions not made for R, do not swap R columns, just rows
  _add_to(columnIndex, columnIndex + 1);  // useless for R?

  return false;
}

template <class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::Pos_index RU_vine_swap<Master_matrix>::_get_death(Index simplexIndex)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    return _matrix()->_death_val(simplexIndex);
  } else {
    if (!_matrix()->reducedMatrixR_.is_zero_column(simplexIndex))
      return _matrix()->reducedMatrixR_.get_column(simplexIndex).get_pivot();

    if constexpr (Master_matrix::Option_list::has_map_column_container) {
      auto it = _matrix()->pivotToColumnIndex_.find(simplexIndex);
      if (it == _matrix()->pivotToColumnIndex_.end()) return Master_matrix::template get_null_value<Pos_index>();
      return it->second;
    } else {
      if (simplexIndex >= _matrix()->pivotToColumnIndex_.size())
        return Master_matrix::template get_null_value<Pos_index>();
      return _matrix()->pivotToColumnIndex_[simplexIndex];
    }
  }
}

template <class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::Pos_index RU_vine_swap<Master_matrix>::_get_birth(
    Index negativeSimplexIndex)
{
  if constexpr (Master_matrix::Option_list::has_column_pairings) {
    return _matrix()->_birth_val(negativeSimplexIndex);
  } else {
    return _matrix()->reducedMatrixR_.get_pivot(negativeSimplexIndex);
  }
}

template <class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::ID_index RU_vine_swap<Master_matrix>::_get_row_id_from_position(
    Pos_index position) const
{
  const auto& map = _matrix()->positionToID_;
  auto it = map.find(position);
  return it == map.end() ? position : it->second;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_RU_VINE_SWAP_H
