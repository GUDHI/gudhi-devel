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

#ifndef PM_BASE_PAIRING_H
#define PM_BASE_PAIRING_H

#include <utility>  //std::swap & std::move
#include <unordered_map>
#include <algorithm>
#include <vector>

#include <gudhi/Persistence_matrix/index_mapper.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Base_pairing, when the computation of the barcode was not enabled or if the pairing
 * is already managed by the vine update classes.
 */
struct Dummy_base_pairing {
  friend void swap([[maybe_unused]] Dummy_base_pairing& d1, [[maybe_unused]] Dummy_base_pairing& d2) noexcept {}
};

/**
 * @class Base_pairing base_pairing.h gudhi/Persistence_matrix/base_pairing.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the barcode for @ref Boundary_matrix if the option was enabled.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Base_pairing
    : protected std::conditional_t<
          Master_matrix::Option_list::has_removable_columns,
          Index_mapper<std::unordered_map<typename Master_matrix::Pos_index, typename Master_matrix::ID_index> >,
          Dummy_index_mapper>
{
 protected:
  using Pos_index = typename Master_matrix::Pos_index;
  using ID_index = typename Master_matrix::ID_index;
  // PIDM = Position to ID Map
  using PIDM = std::conditional_t<Master_matrix::Option_list::has_removable_columns,
                                  Index_mapper<std::unordered_map<Pos_index, ID_index> >,
                                  Dummy_index_mapper>;

 public:
  using Bar = typename Master_matrix::Bar;                           /**< Bar type. */
  using Barcode = typename Master_matrix::Barcode;                   /**< Barcode type. */
  using Column = typename Master_matrix::Column;                     /**< Column type. */
  using Index = typename Master_matrix::Index;                       /**< Container index type. */

  /**
   * @brief Default constructor.
   */
  Base_pairing();

  /**
   * @brief Reduces the matrix stored in @ref Boundary_matrix and computes the corresponding barcode.
   *
   * @warning The barcode will not be recomputed if the matrix is modified later after calling this method
   * for the first time. So call it only once the matrix is finalized. This behaviour could be changed in the future,
   * if the need is mentioned.
   *
   * @return Const reference to the barcode.
   */
  const Barcode& get_current_barcode();

  /**
   * @brief Swap operator.
   */
  friend void swap(Base_pairing& pairing1, Base_pairing& pairing2) noexcept
  {
    if constexpr (Master_matrix::Option_list::has_removable_columns) {
      swap(static_cast<PIDM&>(pairing1), static_cast<PIDM&>(pairing2));
    }
    pairing1.barcode_.swap(pairing2.barcode_);
    pairing1.deathToBar_.swap(pairing2.deathToBar_);
    pairing1.idToPosition_.swap(pairing2.idToPosition_);
    std::swap(pairing1.isReduced_, pairing2.isReduced_);
  }

 protected:
  void _remove_last(Pos_index columnIndex);
  void _insert_id_position(ID_index id, Pos_index pos);
  void _reset();

 private:
  using Dictionary = typename Master_matrix::Bar_dictionary;
  using Base_matrix = typename Master_matrix::Master_boundary_matrix;

  Barcode barcode_;       /**< Bar container. */
  Dictionary deathToBar_; /**< Map from death index to bar index. */
  bool isReduced_;        /**< True if `_reduce()` was called. */
  /**
   * @brief Map from cell ID to cell position. Only stores a pair if ID != position.
   */
  std::unordered_map<ID_index, Pos_index> idToPosition_;  // TODO: test other map types

  void _reduce();
  void _reduce_column_with(Column& toReduce, Column& toAdd);

  // access to inheriting Boundary_matrix class
  constexpr Base_matrix* _matrix() { return static_cast<Base_matrix*>(this); }

  constexpr const Base_matrix* _matrix() const { return static_cast<const Base_matrix*>(this); }
};

template <class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing() : PIDM(), isReduced_(false)
{}

template <class Master_matrix>
inline const typename Base_pairing<Master_matrix>::Barcode& Base_pairing<Master_matrix>::get_current_barcode()
{
  if (!isReduced_) _reduce();
  return barcode_;
}

template <class Master_matrix>
inline void Base_pairing<Master_matrix>::_reduce()
{
  constexpr const Pos_index nullDeath = Master_matrix::template get_null_value<Pos_index>();
  constexpr const Index nullIndex = Master_matrix::template get_null_value<Index>();
  constexpr const ID_index nullPivot = Master_matrix::template get_null_value<ID_index>();

  std::unordered_map<Index, Index> negativeColumns(_matrix()->get_number_of_columns());

  auto dim = _matrix()->get_max_dimension();
  std::vector<std::vector<Index> > columnsByDim(dim + 1);
  for (auto& v : columnsByDim) v.reserve(_matrix()->get_number_of_columns());
  for (unsigned int i = 0; i < _matrix()->get_number_of_columns(); i++) {
    columnsByDim[dim - _matrix()->get_column_dimension(i)].push_back(i);
  }

  for (const auto& cols : columnsByDim) {
    for (Index i : cols) {
      auto& curr = _matrix()->get_column(i);
      if (curr.is_empty()) {
        if (negativeColumns.find(i) == negativeColumns.end()) {
          barcode_.emplace_back(i, nullDeath, dim);
        }
      } else {
        ID_index pivot = curr.get_pivot();
        auto it = idToPosition_.find(pivot);
        Index pivotColumnNumber = it == idToPosition_.end() ? pivot : it->second;
        auto itNeg = negativeColumns.find(pivotColumnNumber);
        Index pivotKiller = itNeg == negativeColumns.end() ? nullIndex : itNeg->second;

        while (pivot != nullPivot && pivotKiller != nullIndex) {
          _reduce_column_with(curr, _matrix()->get_column(pivotKiller));
          pivot = curr.get_pivot();
          it = idToPosition_.find(pivot);
          pivotColumnNumber = it == idToPosition_.end() ? pivot : it->second;
          itNeg = negativeColumns.find(pivotColumnNumber);
          pivotKiller = itNeg == negativeColumns.end() ? nullIndex : itNeg->second;
        }

        if (pivot != nullPivot) {
          negativeColumns.emplace(pivotColumnNumber, i);
          _matrix()->get_column(pivotColumnNumber).clear();
          barcode_.emplace_back(pivotColumnNumber, i, dim - 1);
        } else {
          curr.clear();
          barcode_.emplace_back(i, nullDeath, dim);
        }
      }
    }
    --dim;
  }

  if constexpr (Master_matrix::Option_list::has_removable_columns) {
    // sort barcode by birth such that a removal is trivial
    std::sort(barcode_.begin(), barcode_.end(), [](const Bar& b1, const Bar& b2) { return b1.birth < b2.birth; });
    // map can only be constructed once barcode is sorted
    for (Index i = 0; i < barcode_.size(); ++i) {
      auto d = barcode_[i].death;
      if (d != Master_matrix::template get_null_value<Pos_index>()) {
        deathToBar_.emplace(d, i);
      }
    }
  }

  isReduced_ = true;
}

template <class Master_matrix>
inline void Base_pairing<Master_matrix>::_reduce_column_with(Column& toReduce, Column& toAdd)
{
  if constexpr (Master_matrix::Option_list::is_z2) {
    toReduce += toAdd;
  } else {
    typename Master_matrix::Element coeff = toAdd.get_pivot_value();
    auto& operators = _matrix()->colSettings_->operators;
    coeff = operators.get_inverse(coeff);
    operators.multiply_inplace(coeff, operators.get_characteristic() - toReduce.get_pivot_value());
    toReduce.multiply_source_and_add(toAdd, coeff);
  }
}

template <class Master_matrix>
inline void Base_pairing<Master_matrix>::_remove_last(Pos_index columnIndex)
{
  static_assert(Master_matrix::Option_list::has_removable_columns, "remove_last not available.");

  if (isReduced_) {
    auto it = deathToBar_.find(columnIndex);

    if (it == deathToBar_.end()) {  // birth
      barcode_.pop_back();          // sorted by birth and columnIndex has to be the highest one
    } else {                        // death
      barcode_[it->second].death = Master_matrix::template get_null_value<Pos_index>();
      deathToBar_.erase(it);
    };
  }

  auto it = PIDM::map_.find(columnIndex);
  if (it != PIDM::map_.end()) {
    idToPosition_.erase(it->second);
    PIDM::map_.erase(it);
  }
}

template <class Master_matrix>
inline void Base_pairing<Master_matrix>::_insert_id_position(ID_index id, Pos_index pos)
{
  idToPosition_.emplace(id, pos);
}

template <class Master_matrix>
inline void Base_pairing<Master_matrix>::_reset()
{
  barcode_.clear();
  deathToBar_.clear();
  isReduced_ = false;
  idToPosition_.clear();
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_BASE_PAIRING_H
