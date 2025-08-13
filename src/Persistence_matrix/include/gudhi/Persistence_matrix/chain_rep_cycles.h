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

  /**
   * @brief Default constructor.
   */
  Chain_representative_cycles() = default;

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
inline void Chain_representative_cycles<Master_matrix>::update_representative_cycles()
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

  birthToCycle_.clear();
  birthToCycle_.resize(nberColumns, Master_matrix::template get_null_value<Index>());
  representativeCycles_.clear();

  for (ID_index i = 0; i < nberColumns; i++) {
    auto& col = _matrix()->get_column(_matrix()->get_column_with_pivot(i));
    if (!col.is_paired() || get_position(i) < get_position(_matrix()->get_pivot(col.get_paired_chain_index()))) {
      Cycle cycle;
      if constexpr (is_well_behaved<Master_matrix::Option_list::column_type>::value) {
        cycle.reserve(col.size());
        for (const auto& c : col) {
          if constexpr (Master_matrix::Option_list::is_z2) {
            cycle.push_back(c.get_row_index());
          } else {
            cycle.push_back({c.get_row_index(), c.get_element()});
          }
        }
      } else {
        auto cont = col.get_content();
        for (Index j = 0; j < cont.size(); ++j) {
          if (cont[j] != 0) {
            if constexpr (Master_matrix::Option_list::is_z2) {
              cycle.push_back(j);
            } else {
              cycle.push_back({j, cont[j]});
            }
          }
        }
      }
      representativeCycles_.push_back(std::move(cycle));
      birthToCycle_[get_position(i)] = representativeCycles_.size() - 1;
    }
  }
}

template <class Master_matrix>
inline const std::vector<typename Chain_representative_cycles<Master_matrix>::Cycle>&
Chain_representative_cycles<Master_matrix>::get_representative_cycles()
{
  if (representativeCycles_.empty()) update_representative_cycles();
  return representativeCycles_;
}

template <class Master_matrix>
inline const typename Chain_representative_cycles<Master_matrix>::Cycle&
Chain_representative_cycles<Master_matrix>::get_representative_cycle(const Bar& bar)
{
  if (representativeCycles_.empty()) update_representative_cycles();
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
