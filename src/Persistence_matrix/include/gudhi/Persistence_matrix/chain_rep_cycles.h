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
 * @file chain_rep_cycles.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Chain_representative_cycles class and
 * @ref Gudhi::persistence_matrix::Dummy_chain_representative_cycles structure.
 */

#ifndef PM_CHAIN_REP_CYCLES_H
#define PM_CHAIN_REP_CYCLES_H

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
 * Inherited instead of @ref Chain_representative_cycles, when the computation of the representative cycles
 * were not enabled.
 */
struct Dummy_chain_representative_cycles {
  friend void swap([[maybe_unused]] Dummy_chain_representative_cycles& d1,
                   [[maybe_unused]] Dummy_chain_representative_cycles& d2) noexcept
  {}
};

/**
 * @class Chain_representative_cycles chain_rep_cycles.h gudhi/Persistence_matrix/chain_rep_cycles.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the representative cycles for @ref Chain_matrix if the option was enabled.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Chain_representative_cycles
{
 public:
  using Bar = typename Master_matrix::Bar;                           /**< Bar type. */
  using Cycle = typename Master_matrix::Cycle;                       /**< Cycle type. */
  using Column_container = typename Master_matrix::Column_container; /**< Column container type. */
  using Index = typename Master_matrix::Index;                       /**< @ref MatIdx index type. */
  using ID_index = typename Master_matrix::ID_index;                 /**< @ref IDIdx index type. */
  using Dimension = typename Master_matrix::Dimension;               /**< Dimension type. */

  /**
   * @brief Default constructor.
   */
  Chain_representative_cycles() = default;

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
   * @note For chain matrices with enabled vine swaps, this method will only be more efficient than
   * @ref update_representative_cycles if not called for too many bars.
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
  friend void swap(Chain_representative_cycles& base1, Chain_representative_cycles& base2) noexcept
  {
    base1.representativeCycles_.swap(base2.representativeCycles_);
    base1.birthToCycle_.swap(base2.birthToCycle_);
  }

 protected:
  void _reset();

 private:
  using Master_chain_matrix = typename Master_matrix::Master_chain_matrix;

  std::vector<Cycle> representativeCycles_; /**< Cycle container. */
  std::vector<Index> birthToCycle_;         /**< Map from birth index to cycle index. */

  // access to inheriting Chain_matrix class
  constexpr Master_chain_matrix* _matrix() { return static_cast<Master_chain_matrix*>(this); }

  constexpr const Master_chain_matrix* _matrix() const { return static_cast<const Master_chain_matrix*>(this); }
};

template <class Master_matrix>
inline void Chain_representative_cycles<Master_matrix>::update_representative_cycles(Dimension dim)
{
  Index nberColumns = _matrix()->get_number_of_columns();
  auto get_position = [&](ID_index pivot) {
    if constexpr (Master_matrix::Option_list::has_vine_update) {
      if constexpr (Master_matrix::Option_list::has_map_column_container) {
        return _matrix()->map_.at(pivot);
      } else {
        return _matrix()->map_[pivot];
      }
    } else {
      return pivot;
    }
  };
  Index nullValue = Master_matrix::template get_null_value<Index>();

  birthToCycle_.clear();
  birthToCycle_.resize(nberColumns, nullValue);
  representativeCycles_.clear();

#ifdef GUDHI_USE_TBB
  Index c = 0;
  for (Index i = 0; i < nberColumns; i++) {
    auto& col = _matrix()->get_column(_matrix()->get_column_with_pivot(i));
    if ((dim == Master_matrix::template get_null_value<Dimension>() || _matrix()->get_column_dimension(i) == dim) &&
        (!col.is_paired() || get_position(i) < get_position(_matrix()->get_pivot(col.get_paired_chain_index())))) {
      birthToCycle_[get_position(i)] = c;
      ++c;
    }
  }

  representativeCycles_.resize(c);
  tbb::parallel_for(static_cast<Index>(0), nberColumns, [&](Index i) {
    auto idx = get_position(i);
    if (birthToCycle_[idx] != nullValue) {
      auto& col = _matrix()->get_column(_matrix()->get_column_with_pivot(i));
      representativeCycles_[birthToCycle_[idx]] =
          Master_matrix::build_cycle_from_range(col.get_non_zero_content_range());
    }
  });
#else
  for (ID_index i = 0; i < nberColumns; i++) {
    auto& col = _matrix()->get_column(_matrix()->get_column_with_pivot(i));
    if ((dim == Master_matrix::template get_null_value<Dimension>() || _matrix()->get_column_dimension(i) == dim) &&
        (!col.is_paired() || get_position(i) < get_position(_matrix()->get_pivot(col.get_paired_chain_index())))) {
      representativeCycles_.push_back(Master_matrix::build_cycle_from_range(col.get_non_zero_content_range()));
      birthToCycle_[get_position(i)] = representativeCycles_.size() - 1;
    }
  }
#endif
}

template <class Master_matrix>
inline void Chain_representative_cycles<Master_matrix>::update_representative_cycle(const Bar& bar)
{
  auto nberColumns = _matrix()->get_number_of_columns();
  auto get_position = [&](ID_index pivot) {
    if constexpr (Master_matrix::Option_list::has_vine_update) {
      if constexpr (Master_matrix::Option_list::has_map_column_container) {
        return _matrix()->map_.at(pivot);
      } else {
        return _matrix()->map_[pivot];
      }
    } else {
      return pivot;
    }
  };

  Index nullValue = Master_matrix::template get_null_value<Index>();

  if (birthToCycle_.size() <= bar.birth) {
    birthToCycle_.resize(bar.birth + 1, nullValue);
  }
  if (birthToCycle_[bar.birth] == nullValue) {
    birthToCycle_[bar.birth] = representativeCycles_.size();
    representativeCycles_.resize(representativeCycles_.size() + 1);
  }

  if constexpr (Master_matrix::Option_list::has_vine_update) {
    for (ID_index i = 0; i < nberColumns; i++) {
      if (get_position(i) == bar.birth) {
        auto& col = _matrix()->get_column(_matrix()->get_column_with_pivot(i));
        representativeCycles_[birthToCycle_[bar.birth]] =
            Master_matrix::build_cycle_from_range(col.get_non_zero_content_range());
      }
    }
  } else {
    auto& col = _matrix()->get_column(_matrix()->get_column_with_pivot(bar.birth));
    representativeCycles_[birthToCycle_[bar.birth]] =
        Master_matrix::build_cycle_from_range(col.get_non_zero_content_range());
  }
}

template <class Master_matrix>
inline const std::vector<typename Chain_representative_cycles<Master_matrix>::Cycle>&
Chain_representative_cycles<Master_matrix>::get_representative_cycles()
{
  return representativeCycles_;
}

template <class Master_matrix>
inline const typename Chain_representative_cycles<Master_matrix>::Cycle&
Chain_representative_cycles<Master_matrix>::get_representative_cycle(const Bar& bar)
{
  return representativeCycles_[birthToCycle_[bar.birth]];
}

template <class Master_matrix>
inline void Chain_representative_cycles<Master_matrix>::_reset()
{
  representativeCycles_.clear();
  birthToCycle_.clear();
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_CHAIN_REP_CYCLES_H
