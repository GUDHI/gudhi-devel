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
 * @file ru_rep_cycles.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::RU_representative_cycles class and
 * @ref Gudhi::persistence_matrix::Dummy_ru_representative_cycles structure.
 */

#ifndef PM_RU_REP_CYCLES_H
#define PM_RU_REP_CYCLES_H

#include <utility>    //std::move
#include <algorithm>  //std::sort
#include <vector>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_for.h>
#endif

#include <gudhi/persistence_matrix_options.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref RU_representative_cycles, when the computation of the representative cycles
 * were not enabled.
 */
struct Dummy_ru_representative_cycles {
  friend void swap([[maybe_unused]] Dummy_ru_representative_cycles& d1,
                   [[maybe_unused]] Dummy_ru_representative_cycles& d2) noexcept
  {}
};

/**
 * @class RU_representative_cycles ru_rep_cycles.h gudhi/Persistence_matrix/ru_rep_cycles.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the representative cycles for @ref RU_matrix if the option was enabled.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class RU_representative_cycles
{
 public:
  using Index = typename Master_matrix::Index;         /**< @ref MatIdx index type. */
  using Bar = typename Master_matrix::Bar;             /**< Bar type. */
  using Cycle = typename Master_matrix::Cycle;         /**< Cycle type. */
  using Dimension = typename Master_matrix::Dimension; /**< Dimension type. */

  /**
   * @brief Default constructor.
   */
  RU_representative_cycles() = default;

  /**
   * @brief Computes the current representative cycles of the matrix.
   *
   * @param dim If different from default value, only the cycles of the given dimension are updated.
   * All others are erased.
   */
  void update_representative_cycles(Dimension dim = Master_matrix::template get_null_value<Dimension>());

  /**
   * @brief Computes the current representative cycle of the given bar. All other cycles already computed are left
   * untouched (and therefore they could be unvalid for the current matrix).
   *
   * @param bar Bar corresponding to the wanted representative cycle.
   */
  void update_representative_cycle(const Bar& bar);

  /**
   * @brief Returns the current representative cycles. If the matrix was modified since the last call,
   * @ref update_representative_cycles has to be called to update the returned cycles.
   *
   * @return A const reference to a vector of @ref Matrix::Cycle containing all representative cycles.
   */
  const std::vector<Cycle>& get_representative_cycles();
  /**
   * @brief Returns the representative cycle corresponding to the given bar.
   * If the matrix was modified since the last call, @ref update_representative_cycles or
   * @ref update_representative_cycle has to be called to update the returned cycle.
   *
   * @param bar Bar corresponding to the wanted representative cycle.
   * @return A const reference to the representative cycle.
   */
  const Cycle& get_representative_cycle(const Bar& bar);

  /**
   * @brief Swap operator.
   */
  friend void swap(RU_representative_cycles& base1, RU_representative_cycles& base2) noexcept
  {
    base1.representativeCycles_.swap(base2.representativeCycles_);
    base1.birthToCycle_.swap(base2.birthToCycle_);
  }

 protected:
  void _reset();

 private:
  using Master_RU_matrix = typename Master_matrix::Master_RU_matrix;
  using Inverse_column = Cycle;
  using Content_range = typename Master_matrix::Column::Content_range;

  std::vector<Cycle> representativeCycles_; /**< Cycle container. */
  std::vector<Index> birthToCycle_;         /**< Map from birth index to cycle index. */

  constexpr Master_RU_matrix* _matrix() { return static_cast<Master_RU_matrix*>(this); }
  constexpr const Master_RU_matrix* _matrix() const { return static_cast<const Master_RU_matrix*>(this); }

  void _retrieve_cycle_from_r(Index colIdx, Index repIdx);
  void _retrieve_cycle_from_u(Index colIdx, Index repIdx);
  Inverse_column _get_inverse(Index c);
};

template <class Master_matrix>
inline void RU_representative_cycles<Master_matrix>::update_representative_cycles(Dimension dim)
{
  Index nberColumns = _matrix()->reducedMatrixR_.get_number_of_columns();
  Index nullValue = Master_matrix::template get_null_value<Index>();
  representativeCycles_.clear();
  birthToCycle_.clear();
  birthToCycle_.resize(nberColumns, nullValue);

  Index c = 0;
  for (Index i = 0; i < nberColumns; i++) {
    if ((dim == Master_matrix::template get_null_value<Dimension>() ||
         _matrix()->reducedMatrixR_.get_column_dimension(i) == dim) &&
        _matrix()->reducedMatrixR_.is_zero_column(i)) {
      birthToCycle_[i] = c;
      ++c;
    }
  }

  representativeCycles_.resize(c);
#ifdef GUDHI_USE_TBB
  tbb::parallel_for(static_cast<Index>(0), nberColumns, [&](Index i) {
    if (birthToCycle_[i] != nullValue) {
      Index colIdx = _matrix()->_get_column_with_pivot(i);
      if (colIdx == nullValue) {
        _retrieve_cycle_from_u(i, birthToCycle_[i]);
      } else {
        _retrieve_cycle_from_r(colIdx, birthToCycle_[i]);
      }
    }
  });
#else
  for (Index i = 0; i < nberColumns; ++i) {
    if (birthToCycle_[i] != nullValue) {
      Index colIdx = _matrix()->_get_column_with_pivot(i);
      if (colIdx == nullValue) {
        _retrieve_cycle_from_u(i, birthToCycle_[i]);
      } else {
        _retrieve_cycle_from_r(colIdx, birthToCycle_[i]);
      }
    }
  }
#endif
}

template <class Master_matrix>
inline void RU_representative_cycles<Master_matrix>::update_representative_cycle(const Bar& bar)
{
  Index nullValue = Master_matrix::template get_null_value<Index>();

  if (birthToCycle_.size() <= bar.birth) {
    birthToCycle_.resize(bar.birth + 1, nullValue);
  }
  if (birthToCycle_[bar.birth] == nullValue) {
    birthToCycle_[bar.birth] = representativeCycles_.size();
    representativeCycles_.resize(representativeCycles_.size() + 1);
  }

  Index colIdx = _matrix()->_get_column_with_pivot(bar.birth);
  if (colIdx == nullValue) {
    _retrieve_cycle_from_u(bar.birth, birthToCycle_[bar.birth]);
  } else {
    _retrieve_cycle_from_r(colIdx, birthToCycle_[bar.birth]);
  }
}

template <class Master_matrix>
inline const std::vector<typename RU_representative_cycles<Master_matrix>::Cycle>&
RU_representative_cycles<Master_matrix>::get_representative_cycles()
{
  return representativeCycles_;
}

template <class Master_matrix>
inline const typename RU_representative_cycles<Master_matrix>::Cycle&
RU_representative_cycles<Master_matrix>::get_representative_cycle(const Bar& bar)
{
  return representativeCycles_[birthToCycle_[bar.birth]];
}

template <class Master_matrix>
inline void RU_representative_cycles<Master_matrix>::_retrieve_cycle_from_r(Index colIdx, Index repIdx)
{
  auto& col = _matrix()->reducedMatrixR_.get_column(colIdx);
  representativeCycles_[repIdx] = Master_matrix::build_cycle_from_range(col.get_non_zero_content_range());
}

template <class Master_matrix>
inline void RU_representative_cycles<Master_matrix>::_retrieve_cycle_from_u(Index colIdx, Index repIdx)
{
  // TODO: if rep_cycles true but not vineyards, this could be avoided by directly computing V instead of U
  representativeCycles_[repIdx] = _get_inverse(colIdx);
}

template <class Master_matrix>
inline typename RU_representative_cycles<Master_matrix>::Inverse_column
RU_representative_cycles<Master_matrix>::_get_inverse(Index c)
{
  using E = typename Master_matrix::Element;
  auto& matrix = _matrix()->mirrorMatrixU_;
  auto size = matrix.get_number_of_columns();
  [[maybe_unused]] const auto& op = _matrix()->operators_;
  Inverse_column res;

  auto _last_diagonal_value = [&]() -> E {
    auto& col = matrix.get_column(size - 1);
    if (col.is_empty()) return 0;  // happens only if the user multiplied by 0 a column in R
    // if column not empty, pivot has to be on the diagonal
    if constexpr (Master_matrix::Option_list::is_z2) {
      return 1;
    } else {
      return op->get_inverse(matrix.get_column(size - 1).get_pivot_value());
    }
  };

  auto _push_cell = [&](auto i, E e) -> void {
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (e) res.push_back(i);
    } else {
      if (e != op->get_additive_identity()) res.push_back({i, e});
    }
  };

  auto _substract = [&](E& e, auto resIt, const auto& cell) -> void {
    if (resIt != res.rend() && Master_matrix::get_row_index(*resIt) == cell.get_row_index()) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        e = !e;
      } else {
        op->subtract_inplace_front(e, cell.get_element() * Master_matrix::get_element(*resIt));
      }
    }
  };

  auto _multiply = [&](E& e, E m) -> void {
    if constexpr (Master_matrix::Option_list::is_z2) {
      e = m && e;  // just in case, but m is only possibly 0 if the user multiplied by 0 a column in R
    } else {
      op->multiply_inplace(e, op->get_inverse(m));
    }
  };

  auto _translate = [&](std::size_t i) -> void {
    const auto& map = _matrix()->positionToID_;
    auto& idx = Master_matrix::get_row_index(res[i]);
    auto it = map.find(idx);
    if (it != map.end()) idx = it->second;
  };

  if (c == size - 1) _push_cell(size - 1, _last_diagonal_value());
  for (int i = size - 2; i >= 0; --i) {
    E e = static_cast<int>(c) == i;
    auto& line = matrix.get_column(i);
    Content_range r = line.get_non_zero_content_range();
    auto resIt = res.rbegin();
    auto lineIt = r.begin();
    E diag(0);
    if (static_cast<int>(lineIt->get_row_index()) == i) {
      diag = lineIt->get_element();
      ++lineIt;
    }
    while (lineIt != r.end() && resIt != res.rend()) {
      while (resIt != res.rend() && Master_matrix::get_row_index(*resIt) < lineIt->get_row_index()) ++resIt;
      _substract(e, resIt, *lineIt);
      ++lineIt;
    }
    _multiply(e, diag);
    _push_cell(i, e);
  }

  // reverse order + PosIdx to IDIdx translation
  for (std::size_t incr = 0, decr = res.size() - 1; incr < decr; ++incr, --decr) {
    _translate(incr);
    _translate(decr);
    std::swap(res[incr], res[decr]);
  }

  return res;
}

template <class Master_matrix>
inline void RU_representative_cycles<Master_matrix>::_reset()
{
  representativeCycles_.clear();
  birthToCycle_.clear();
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_RU_REP_CYCLES_H
