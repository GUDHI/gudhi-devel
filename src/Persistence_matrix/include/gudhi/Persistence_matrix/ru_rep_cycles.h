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
struct Dummy_ru_representative_cycles
{
  friend void swap([[maybe_unused]] Dummy_ru_representative_cycles& d1,
                   [[maybe_unused]] Dummy_ru_representative_cycles& d2) {}
};

// TODO: add coefficients ? Only Z2 token into account for now.
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
  using Index = typename Master_matrix::Index;  /**< @ref MatIdx index type. */
  using Bar = typename Master_matrix::Bar;      /**< Bar type. */
  using Cycle = typename Master_matrix::Cycle;  /**< Cycle type. */

  /**
   * @brief Default constructor.
   */
  RU_representative_cycles();
  /**
   * @brief Copy constructor.
   * 
   * @param matrixToCopy Matrix to copy.
   */
  RU_representative_cycles(const RU_representative_cycles& matrixToCopy);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  RU_representative_cycles(RU_representative_cycles&& other) noexcept;

  /**
   * @brief Computes the current representative cycles of the matrix.
   */
  void update_representative_cycles();

  /**
   * @brief Returns the current representative cycles. If the matrix is modified later after the first call,
   * @ref update_representative_cycles has to be called to update the returned cycles.
   * 
   * @return A const reference to a vector of @ref Matrix::Cycle containing all representative cycles.
   */
  const std::vector<Cycle>& get_representative_cycles();
  /**
   * @brief Returns the representative cycle corresponding to the given bar.
   * If the matrix is modified later after the first call,
   * @ref update_representative_cycles has to be called to update the returned cycles.
   * 
   * @param bar Bar corresponding to the wanted representative cycle.
   * @return A const reference to the representative cycle.
   */
  const Cycle& get_representative_cycle(const Bar& bar);

  /**
   * @brief Assign operator.
   */
  RU_representative_cycles& operator=(RU_representative_cycles other);
  /**
   * @brief Swap operator.
   */
  friend void swap(RU_representative_cycles& base1, RU_representative_cycles& base2) {
    base1.representativeCycles_.swap(base2.representativeCycles_);
    base1.birthToCycle_.swap(base2.birthToCycle_);
  }

 private:
  using Master_RU_matrix = typename Master_matrix::Master_RU_matrix;

  std::vector<Cycle> representativeCycles_; /**< Cycle container. */
  std::vector<Index> birthToCycle_;         /**< Map from birth index to cycle index. */

  constexpr Master_RU_matrix* _matrix() { return static_cast<Master_RU_matrix*>(this); }
  constexpr const Master_RU_matrix* _matrix() const { return static_cast<const Master_RU_matrix*>(this); }

  // TODO: if rep_cycles true but not vineyards, this could be avoided by directly computing V instead of U
  std::vector<std::vector<typename Master_matrix::Entry_representative>> _get_inverse();
};

template <class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles() 
{}

template <class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles(
    const RU_representative_cycles<Master_matrix>& matrixToCopy)
    : representativeCycles_(matrixToCopy.representativeCycles_), birthToCycle_(matrixToCopy.birthToCycle_) 
{}

template <class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles(
    RU_representative_cycles<Master_matrix>&& other) noexcept
    : representativeCycles_(std::move(other.representativeCycles_)), birthToCycle_(std::move(other.birthToCycle_)) 
{}

template <class Master_matrix>
inline void RU_representative_cycles<Master_matrix>::update_representative_cycles() 
{
  birthToCycle_.clear();
  birthToCycle_.resize(_matrix()->reducedMatrixR_.get_number_of_columns(), -1);
  Index c = 0;
  for (Index i = 0; i < _matrix()->reducedMatrixR_.get_number_of_columns(); i++) {
    if (_matrix()->reducedMatrixR_.is_zero_column(i)) {
      birthToCycle_[i] = c;
      ++c;
    }
  }
  representativeCycles_.clear();
  representativeCycles_.resize(c);
  for (Index i = 0; i < _matrix()->mirrorMatrixU_.get_number_of_columns(); i++) {
    for (const auto& entry : _matrix()->mirrorMatrixU_.get_column(i)) {
      auto idx = birthToCycle_[entry.get_row_index()];
      if (idx != Master_matrix::template get_null_value<Index>()) {
        representativeCycles_[idx].push_back(i);
      }
    }
  }
}

template <class Master_matrix>
inline const std::vector<typename RU_representative_cycles<Master_matrix>::Cycle>&
RU_representative_cycles<Master_matrix>::get_representative_cycles() 
{
  if (representativeCycles_.empty()) update_representative_cycles();
  return representativeCycles_;
}

template <class Master_matrix>
inline const typename RU_representative_cycles<Master_matrix>::Cycle&
RU_representative_cycles<Master_matrix>::get_representative_cycle(const Bar& bar) 
{
  if (representativeCycles_.empty()) update_representative_cycles();
  return representativeCycles_[birthToCycle_[bar.birth]];
}

template <class Master_matrix>
inline RU_representative_cycles<Master_matrix>& RU_representative_cycles<Master_matrix>::operator=(
    RU_representative_cycles<Master_matrix> other) 
{
  representativeCycles_.swap(other.representativeCycles_);
  birthToCycle_.swap(other.birthToCycle_);
  return *this;
}

template <class Master_matrix>
inline std::vector<std::vector<typename Master_matrix::Entry_representative> >
RU_representative_cycles<Master_matrix>::_get_inverse()
{
  using E = typename Master_matrix::Element;
  auto& matrix = _matrix()->mirrorMatrixU_;
  auto size = matrix.get_number_of_columns();
  [[maybe_unused]] const auto& op = _matrix()->operators_;
  std::vector<std::vector<typename Master_matrix::Entry_representative>> res(size);

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

  auto _push_cell = [&](int c, int i, E e) -> void {
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (e) res[c].push_back(i);
    } else {
      if (e != op->get_additive_identity()) res[c].push_back({i, e});
    }
  };

  auto _substract = [&](E& e, auto resIt, const auto& cell, int c) -> void {
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (resIt != res[c].rend() && *resIt == cell.get_row_index()) e = !e;
    } else {
      if (resIt != res[c].rend() && resIt->first == cell.get_row_index())
        op->subtract_inplace_front(e, cell.get_element() * resIt->second);
    }
  };

  auto _substract_vec = [&](E& e, auto p, auto line) -> void {
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (line[p]) e = !e;
    } else {
      if (line[p.first]) op->subtract_inplace_front(e, line[p.first] * p.second);
    }
  };

  auto _multiply = [&](E& e, E m) -> void {
    if constexpr (Master_matrix::Option_list::is_z2) {
      e = m && e;  // just in case, but m is only possibly 0 if the user multiplied by 0 a column in R
    } else {
      op->multiply_inplace(e, op->get_inverse(m));
    }
  };

  auto _assign = [&](E& e, const auto& cell) -> void {
    if constexpr (Master_matrix::Option_list::is_z2) {
      e = !e;
    } else {
      e = cell.get_element();
    }
  };

  auto _get_index = [&](auto resIt) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      return *resIt;
    } else {
      return resIt->first;
    }
  };

  _push_cell(size - 1, size - 1, _last_diagonal_value());
  for (int c = size - 1; c >= 0; --c) {
    for (int i = size - 2; i >= 0; --i) {
      E e = c == i ? 1 : 0;
      // ugly......
      if constexpr (is_well_behaved<Master_matrix::Option_list::column_type>::value) {
        const auto& line = matrix.get_column(i);
        auto resIt = res[c].rbegin();
        auto lineIt = line.begin();
        E diag(0);
        if (static_cast<int>(lineIt->get_row_index()) == i) {
          _assign(diag, *lineIt);
          ++lineIt;
        }
        while (lineIt != line.end() && resIt != res[c].rend()) {
          while (resIt != res[c].rend() && _get_index(resIt) < lineIt->get_row_index()) ++resIt;
          _substract(e, resIt, *lineIt, c);
          ++lineIt;
        }
        _multiply(e, diag);
      } else {
        auto line = matrix.get_column(i).get_content(size); // linear...
        for (const auto& p : res[c]) {
          _substract_vec(e, p, line);
        }
        _multiply(e, line[i]);
      }
      _push_cell(c, i, e);
    }
    std::reverse(res[c].begin(), res[c].end());
  }

  return res;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_RU_REP_CYCLES_H
