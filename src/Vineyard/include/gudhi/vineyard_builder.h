/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file vineyard_builder.h
 * @author Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::vineyard::Vineyard_builder class.
 */

#ifndef GUDHI_VINEYARD_BUILDER_H_
#define GUDHI_VINEYARD_BUILDER_H_

#include <array>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <optional>

#include <gudhi/Debug_utils.h>
#include <gudhi/simple_mdspan.h>
#include <gudhi/vineyard_base.h>
#include <gudhi/persistence_interval.h>

namespace Gudhi {
namespace vineyard {

template <typename T, typename Dimension>
class Vine
{
 public:
  using Coordinate = std::array<T, 2>;
  using Coordinate_container = std::vector<Coordinate>;

  Vine(Dimension dim, int numberOfUpdates = 1) : dim_(dim) { coordinates_.reserve(numberOfUpdates); }

  void add_pair(T birth, T death) { coordinates_.push_back({birth, death}); }

  template <typename Index>
  const Coordinate& get_pair(Index i) const
  {
    return coordinates_[i];
  }

  const Coordinate_container& get_pairs() const { return coordinates_; }

 private:
  Dimension dim_;
  Coordinate_container coordinates_;
};

/**
 * @class Vineyard_builder vineyard_builder.h gudhi/vineyard_builder.h
 * @brief
 * @details
 *
 * @ingroup vineyard
 *
 * @tparam VineyardOptions Structure following the @ref VineyardOptions concept. Default value: @ref
 * Default_vineyard_options.
 */
template <typename T, class VineyardOptions = Default_vineyard_options, bool flat = false>
class Vineyard_builder
{
 public:
  using value_type = T;
  using Base = Vineyard_base<VineyardOptions>;
  using Index = typename Base::Index;                                        /**< Complex index type. */
  using Dimension = typename Base::Dimension;                                /**< Dimension type. */
  using Bar = Gudhi::persistence_matrix::Persistence_interval<Dimension, T>; /**< Bar type. */
  using Cycle = typename Base::Cycle;                                        /**< Cycle type. */
  using Vine_t = Vine<T, Dimension>;
  using Flat_vines = std::vector<std::array<T, 2> >;
  using Vineyard = std::conditional_t<flat, std::vector<Flat_vines>, std::vector<Vine_t> >;

  Vineyard_builder(bool storeRepCycles = false, Dimension repCyclesDim = Base::nullDimension)
  {
    if (storeRepCycles) {
      latest_representative_cycles_.emplace();
      repCyclesDim_.emplace(repCyclesDim);
    }
  }

  template <class Boundary_range, class Dimension_range, class Filtration_range>
  void initialize(const Boundary_range& boundaryMatrix,
                  const Dimension_range& dimensions,
                  const Filtration_range& filtrationValues,
                  int numberOfUpdates = 0)
  {
    base_.initialize(boundaryMatrix, dimensions, filtrationValues);
    const auto& barcode = base_.get_current_barcode();  // forward only range
    vineyard_.clear();
    numberOfBars_.clear();
    if (latest_representative_cycles_) {
      latest_representative_cycles_->clear();
    }
    Index idx = 0;
    if constexpr (flat) {
      std::vector<Index> latestIndex;
      for (const auto& bar : barcode) {
        if (bar.dim >= static_cast<Dimension>(vineyard_.size())) {
          vineyard_.resize(bar.dim + 1);
          latestIndex.resize(bar.dim + 1, -1);
        }
        ++latestIndex[bar.dim];
        vineyard_[bar.dim].push_back(
            {filtrationValues[bar.birth], bar.death == Base::Bar::inf ? Bar::inf : filtrationValues[bar.death]});
        if (_store_cycle(bar, filtrationValues)) {
          auto cycle = base_.get_current_representative_cycle(idx, true);
          latest_representative_cycles_->emplace_back(Cycle(cycle.begin(), cycle.end()), bar.dim, latestIndex[bar.dim]);
        }
        ++idx;
      }
      numberOfBars_.resize(vineyard_.size());
      for (Index i = 0; i < vineyard_.size(); ++i) {
        numberOfBars_[i] = vineyard_[i].size();
        vineyard_[i].reserve(numberOfBars_[i] * (numberOfUpdates + 1));
      }
    } else {
      for (const auto& bar : barcode) {
        vineyard_.emplace_back(bar.dim, numberOfUpdates + 1);
        vineyard_.back().add_pair(filtrationValues[bar.birth],
                                  bar.death == Base::Bar::inf ? Bar::inf : filtrationValues[bar.death]);
        if (bar.dim >= numberOfBars_.size()) numberOfBars_.resize(bar.dim + 1, 0);
        ++numberOfBars_[bar.dim];
        if (_store_cycle(bar, filtrationValues)) {
          auto cycle = base_.get_current_representative_cycle(idx, true);
          latest_representative_cycles_->emplace_back(Cycle(cycle.begin(), cycle.end()), bar.dim, idx);
        }
        ++idx;
      }
    }
  }

  template <class Filtration_range>
  void update(const Filtration_range& filtrationValues)
  {
    base_.update(filtrationValues);
    const auto& barcode = base_.get_current_barcode();  // forward only range + order is preserved
    Index idx = 0;
    std::vector<Index> latestIndex;
    if constexpr (flat) latestIndex.resize(vineyard_.size(), -1);
    if (latest_representative_cycles_) {
      latest_representative_cycles_->clear();
    }
    for (const auto& bar : barcode) {
      T birth = filtrationValues[bar.birth];
      T death = bar.death == Base::Bar::inf ? Bar::inf : filtrationValues[bar.death];
      if constexpr (flat) {
        ++latestIndex[bar.dim];
        vineyard_[bar.dim].push_back({birth, death});
      } else {
        vineyard_[idx].add_pair(birth, death);
      }
      if (_store_cycle(bar, filtrationValues)) {
        const auto& cycle = base_.get_current_representative_cycle(idx, true);
        latest_representative_cycles_->emplace_back(
            Cycle(cycle.begin(), cycle.end()), bar.dim, flat ? latestIndex[bar.dim] : idx);
      }
      ++idx;
    }
  }

  const Vineyard& get_current_vineyard() const { return vineyard_; }

  const std::vector<Index>& get_number_of_vines_by_dimension() const { return numberOfBars_; }

  const std::vector<std::tuple<Cycle, Dimension, Index> >& get_latest_representative_cycles()
  {
    if (latest_representative_cycles_) return *latest_representative_cycles_;
    throw std::invalid_argument("Representative cycles were not stored.");
  }

 private:
  Base base_;
  Vineyard vineyard_;
  std::vector<Index> numberOfBars_;
  std::optional<Dimension> repCyclesDim_;
  std::optional<std::vector<std::tuple<Cycle, Dimension, Index> > > latest_representative_cycles_;

  template <class Filtration_range>
  bool _store_cycle(const typename Base::Bar& bar, const Filtration_range& filtrationValues) const
  {
    if (!repCyclesDim_) return false;
    if (*repCyclesDim_ != Base::nullDimension && bar.dim != *repCyclesDim_) return false;
    if (filtrationValues[bar.birth] == Bar::inf) return false;
    return bar.death == Base::Bar::inf || filtrationValues[bar.death] - filtrationValues[bar.birth] > 0;
  }
};

}  // namespace vineyard
}  // namespace Gudhi

#endif  // GUDHI_VINEYARD_BUILDER_H_
