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
 * @brief Contains the implementation of the @ref Gudhi::vineyard::Vineyard_builder class, as well as the
 * @ref Gudhi::vineyard::Vine class.
 */

#ifndef GUDHI_VINEYARD_BUILDER_H_
#define GUDHI_VINEYARD_BUILDER_H_

#include <array>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <optional>

#include <gudhi/Debug_utils.h>
#include <gudhi/vineyard_base.h>
#include <gudhi/persistence_interval.h>

namespace Gudhi {
namespace vineyard {

/**
 * @ingroup vineyard
 *
 * @brief Class describing a vine.
 *
 * @tparam T Birth and death value type.
 * @tparam Dimension Dimension type. Has to be an integer type.
 */
template <typename T, typename Dimension>
class Vine {
 public:
  using Coordinate = std::array<T, 2>;                  /**< Vine coordinate type. Array of 2 elements. */
  using Coordinate_container = std::vector<Coordinate>; /**< Vector of Coordinate. */

  /**
   * @brief Constructs an vine of given dimension.
   *
   * @param dim Dimension of the vine.
   * @param numberOfUpdates Completely optional. Predicted length of the vine to preallocate memory for it.
   * Default value: 1.
   */
  Vine(Dimension dim, int numberOfUpdates = 1) : dim_(dim) { coordinates_.reserve(numberOfUpdates); }

  /**
   * @brief Adds a birth/death pair to the end of the vine, i.e., its next coordinates.
   */
  void add_pair(T birth, T death) { coordinates_.push_back({birth, death}); }

  /**
   * @brief Returns the birth/death pair at given index in the vine.
   *
   * @tparam Index Integer type.
   */
  template <typename Index>
  const Coordinate& get_pair(Index i) const {
    return coordinates_[i];
  }

  /**
   * @brief Returns the death/birth pair container.
   */
  const Coordinate_container& get_pairs() const { return coordinates_; }

  /**
   * @brief Returns the dimension of the vine.
   */
  Dimension get_dimension() const { return dim_; }

  /**
   * @brief Returns the number of pairs in the vine.
   */
  std::size_t size() const { return coordinates_.size(); }

  /**
   * @brief Basic equality operator
   */
  friend bool operator==(const Vine& v1, const Vine& v2) {
    return v1.dim_ == v2.dim_ && v1.coordinates_ == v2.coordinates_;
  }

  /**
   * @brief Basic unequality operator
   */
  friend bool operator!=(const Vine& v1, const Vine& v2) { return !(v1 == v2); }

  friend std::ostream& operator<<(std::ostream& stream, const Vine& vine)
  {
    stream << "[" << vine.dim_ << "] ";
    for (const auto& v : vine.coordinates_) {
      if constexpr (std::numeric_limits<T>::has_infinity) {
        stream << "(" << v[0] << ", " << v[1] << ") ";
      } else {
        stream << "(";
        if (v[0] == std::numeric_limits<T>::max())
          stream << "inf";
        else
          stream << v[0];
        stream << " - ";
        if (v[1] == std::numeric_limits<T>::max())
          stream << "inf";
        else
          stream << v[1];
        stream << ") ";
      }
    }
    return stream;
  }

 private:
  Dimension dim_;
  Coordinate_container coordinates_;
};

/**
 * @class Vineyard_builder vineyard_builder.h gudhi/vineyard_builder.h
 *
 * @brief Class building and storing a vineyard from given filtered complex and filtration updates. If desired, also
 * stores the representative cycles from the latest step (all or just particular dimension).
 *
 * @details After constructing the class with the desired options, the vineyard is build by first calling
 * @ref Vineyard_builder::initialize to build the first barcode. The remaining of the vines are then computed by
 * calling @ref Vineyard_builder::update in the right order with the new filtration values. Once all updates done,
 * the completed vineyard can be retrieved with @ref Vineyard_builder::get_current_vineyard.
 *
 * @ingroup vineyard
 *
 * @tparam T Type of the filtration values.
 * @tparam VineyardOptions Structure following the @ref VineyardOptions concept. Default value:
 * @ref Default_vineyard_options.
 * @tparam flat Determines the %Vine and Vineyard type. See @ref Vineyard_builder::Vineyard for more information.
 */
template <typename T, class VineyardOptions = Default_vineyard_options, bool flat = false>
class Vineyard_builder {
 public:
  using value_type = T;                                                      /**< Filtration value type. */
  using Base = Vineyard_base<VineyardOptions>;                               /**< Base computing the actual updates. */
  using Index = typename Base::Index;                                        /**< Complex index type. */
  using Dimension = typename Base::Dimension;                                /**< Dimension type. */
  using Bar = Gudhi::persistence_matrix::Persistence_interval<Dimension, T>; /**< Bar type. */
  using Cycle = typename Base::Cycle;                                        /**< Cycle type. */
  /**
   * @brief %Vine type if the template `flat` of the class is set to false: @ref Vine with template arguments deduced
   * from template arguments `T` and `VineyardOptions`.
   */
  using Vine_t = Vine<T, Dimension>;
  /**
   * @brief Type used for @ref Vineyard if the template `flat` of the class is set to true. Representation for all
   * vines of same dimension in form of a "tensor": vector of 2-arrays (i.e. contiguous in memory) such that the
   * \f$ n \f$ first element represent the first coordinates of all \f$ n \f$ vines in the vineyard, the \f$ n \f$
   * next elements the second coordinates of the vines and so on... I.e. it is a flat tensor of shape
   * \f$ (update\_number \times vine\_number \times 2) \f$. The number \f$ n \f$ of vines can be retrieved with the
   * method @ref get_number_of_vines_by_dimension.
   */
  using Flat_vines = std::vector<std::array<T, 2> >;  // used for python
  /**
   * @brief Vineyard type. If the template `flat` of the class is set to false, it is a vector of @ref Vine_t
   * (i.e. @ref Vine). If the template is set to true, it is a vector of @ref Flat_vines where the vines at index
   * \f$ d \f$ are the vines of dimension \f$ d \f$.
   */
  using Vineyard = std::conditional_t<flat, std::vector<Flat_vines>, std::vector<Vine_t> >;

  /**
   * @brief Constructor.
   *
   * @param storeRepCycles If true, @ref get_latest_representative_cycles will not throw and representative cycles will
   * be stored depending on the next argument `repCyclesDim`.
   * @param repCyclesDim Ignored if @p storeRepCycles is false. If set with a positive integer, only representative
   * cycles of this dimension are stored. If unset (i.e. set with default value), all cycles are stored. If set with
   * any other value, nothing is stored.
   */
  Vineyard_builder(bool storeRepCycles = false, Dimension repCyclesDim = Base::nullDimension) {
    if (storeRepCycles) {
      latest_representative_cycles_.emplace();
      repCyclesDim_.emplace(repCyclesDim);
    }
  }

  /**
   * @brief Initializes the vineyard with the first barcode recomputed from scratch. Any vineyard or representative
   * cycle computed before will be cleared and replaced.
   *
   * @tparam BoundaryRange Range of ranges of integers. Has to implement a `size` method and an `operator[]` with
   * a another nested `operator[]`.
   * @tparam DimensionRange Range of integers. Has to implement a `operator[]` method.
   * @tparam FiltrationRange Range of arithmetic values or at least of values with an `operator<`. Has to implement a
   * `operator[]` method.
   * @param boundaryMatrix Boundary container of the filtered complex. The container does not need to be ordered, but
   * the boundaries have to be represented by the indices of their faces in the container.
   * @param dimensions Dimension container of cells in the complex. A value at index \f$ i \f$ has to correspond to the
   * dimension of the cell represented by the boundary at index \f$ i \f$ in `boundaryMatrix`.
   * @param filtrationValues Filtration value container of cells in the complex. A value at index \f$ i \f$ has to
   * correspond to the filtration value of the cell represented by the boundary at index \f$ i \f$ in `boundaryMatrix`.
   * Note that the filtration is assumed to be a 1-parameter filtration.
   * @param numberOfUpdates Optional. Predicted number of updates after this initialization, to preallocate memory for
   * the vines. Default value: 0.
   */
  template <class BoundaryRange, class DimensionRange, class FiltrationRange,
            class = std::enable_if_t<Gudhi::persistence_matrix::RangeTraits<FiltrationRange>::has_size> >
  void initialize(const BoundaryRange& boundaryMatrix, const DimensionRange& dimensions,
                  const FiltrationRange& filtrationValues, int numberOfUpdates = 0) {
    base_.initialize(boundaryMatrix, dimensions, filtrationValues);
    _initialize_stored_vineyard(
        numberOfUpdates, [&](auto b) -> T { return filtrationValues[b]; },
        [&](auto d) -> T { return d == Base::Bar::inf ? Bar::inf : filtrationValues[d]; },
        [&](const auto& bar) -> bool { return _store_cycle(bar, filtrationValues); });
  }

  /**
   * @brief Initializes the vineyard with the first barcode recomputed from scratch. Any vineyard or representative
   * cycle computed before will be cleared and replaced. Different from the other @ref initialize, no filtration value
   * is associated to the cells, so the stored vineyard will use the indices of the faces in the given boundary matrix
   * instead.
   * 
   * @tparam BoundaryRange Range of ranges of integers. Has to implement a `size` method and an `operator[]` with
   * a another nested `operator[]`.
   * @tparam DimensionRange Range of integers. Has to implement a `operator[]` method.
   * @tparam F Method taking to integer values and returning a boolean.
   * @param boundaryMatrix Boundary container of the filtered complex. The container does not need to be ordered, but
   * the boundaries have to be represented by the indices of their faces in the container.
   * @param dimensions Dimension container of cells in the complex. A value at index \f$ i \f$ has to correspond to the
   * dimension of the cell represented by the boundary at index \f$ i \f$ in `boundaryMatrix`.
   * @param is_before_in_filtration Method taking two indices \f$ i \f$ and \f$ j \f$, returning `true` if and only if
   * the cell at index \f$ i \f$ in the boundary matrix should appear in the filtration strictly before the cell
   * at index \f$ j \f$. Note that the stored filtration will first be sorted by dimension, so there is no particular
   * need to distinguish two cells with same filtration value: the order induced by this parameter has to be total
   * for filtration values, but not for the cells them-selves.
   * @param numberOfUpdates Optional. Predicted number of updates after this initialization, to preallocate memory for
   * the vines. Default value: 0.
   */
  template <class BoundaryRange, class DimensionRange, class F,
            class = std::enable_if_t<!Gudhi::persistence_matrix::RangeTraits<F>::has_size> >
  void initialize(const BoundaryRange& boundaryMatrix, const DimensionRange& dimensions, F&& is_before_in_filtration,
                  int numberOfUpdates = 0) {
    base_.initialize(boundaryMatrix, dimensions, std::forward<F>(is_before_in_filtration));
    _initialize_stored_vineyard(
        numberOfUpdates, [&](auto b) -> T { return b; },
        [&](auto d) -> T { return d == Base::Bar::inf ? Bar::inf : static_cast<T>(d); },
        [&](const auto& bar) -> bool { return _store_cycle(bar); });
  }

  /**
   * @brief Extends the vines with the barcode from the new given filtration values. The underlying complex is assumed
   * to be the same than at initialization. Is meant to be used with
   * @ref initialize(const BoundaryRange&,const DimensionRange&,const FiltrationRange&,int).
   *
   * @pre The first barcode has to have been initialized with @ref initialize.
   *
   * @tparam FiltrationRange Range of arithmetic values or at least of values with an `operator<`. Has to implement a
   * `operator[]` method.
   * @param filtrationValues New filtration value container. As at initialization, a value at index \f$ i \f$ has to
   * correspond to the filtration value of the cell represented by the boundary at index \f$ i \f$ in the initializing
   * argument `boundaryMatrix`. Note that the filtration is assumed to be a 1-parameter filtration.
   */
  template <class FiltrationRange, class = std::enable_if_t<!std::is_integral_v<FiltrationRange> > >
  void update(const FiltrationRange& filtrationValues) {
    base_.update(filtrationValues);
    _update_stored_vineyard([&](auto b) -> T { return filtrationValues[b]; },
                            [&](auto d) -> T { return d == Base::Bar::inf ? Bar::inf : filtrationValues[d]; },
                            [&](const auto& bar) -> bool { return _store_cycle(bar, filtrationValues); });
  }

  /**
   * @brief Extends the vines with the barcode generated by the swap of the cell of given filtration index and the
   * cell coming just after in the filtration. The new filtration after the swap has to be a valid filtration,
   * otherwise the behaviour is undefined. Different from @ref update(const FiltrationRange&), the barcode is
   * represented by the indices of the cells in the complex at initialization, not their filtration values.
   * Is meant to be used with @ref initialize(const BoundaryRange&,const DimensionRange&,F&&,int).
   * 
   * @tparam I Integer type.
   * @param i Position in the filtration (not in the complex) of the cell to swap with the cell just after it in
   * the filtration.
   * @param store If false, the new obtained barcode is not stored in the vineyard, only the current state is updated.
   */
  template <typename I, class = std::enable_if_t<std::is_integral_v<I> > >
  void update(I i, bool store = false) {
    base_.update(i);
    if (store) {
      _update_stored_vineyard([&](auto b) -> T { return b; },
                              [&](auto d) -> T { return d == Base::Bar::inf ? Bar::inf : static_cast<T>(d); },
                              [&](const auto& bar) -> bool { return _store_cycle(bar); });
    }
  }

  /**
   * @brief Returns the vineyard in its current state.
   */
  const Vineyard& get_current_vineyard() const { return vineyard_; }

  /**
   * @brief Returns the number of vines in each dimension (dimension 0 at index 0 etc.).
   */
  const std::vector<Index>& get_number_of_vines_by_dimension() const { return numberOfBars_; }

  std::size_t get_filtration_size() const {
    return base_.get_current_order().size();
  }

  /**
   * @brief Returns the representative cycles of the barcode of the last update.
   *
   * A cycle is represented as a tuple of three values: the cycle it-self (range of indices of the cells), the
   * dimension of the cycle and the bar/vine index the cycle corresponds to. In this order.
   *
   * @exception std::invalid_argument If the stored option was set to false at construction.
   */
  const std::vector<std::tuple<Cycle, Dimension, Index> >& get_latest_representative_cycles() const {
    if (latest_representative_cycles_) return *latest_representative_cycles_;
    throw std::invalid_argument("Representative cycles were not stored.");
  }

  static constexpr bool has_flat_vineyard() { return flat; }

 private:
  Base base_;                             /**< Vine computation base. */
  Vineyard vineyard_;                     /**< Vineyard. */
  std::vector<Index> numberOfBars_;       /**< Number fo vines by dimension. */
  std::optional<Dimension> repCyclesDim_; /**< Dimension of stored cycles. */
  std::optional<std::vector<std::tuple<Cycle, Dimension, Index> > > latest_representative_cycles_; /**< Cycles. */

  template <class FiltrationRange>
  bool _store_cycle(const typename Base::Bar& bar, const FiltrationRange& filtrationValues) const {
    if (!repCyclesDim_) return false;
    if (*repCyclesDim_ != Base::nullDimension && bar.dim != *repCyclesDim_) return false;
    if (filtrationValues[bar.birth] == Bar::inf) return false;
    return bar.death == Base::Bar::inf || filtrationValues[bar.death] - filtrationValues[bar.birth] > 0;
  }

  bool _store_cycle(const typename Base::Bar& bar) const {
    if (!repCyclesDim_) return false;
    if (*repCyclesDim_ != Base::nullDimension && bar.dim != *repCyclesDim_) return false;
    return true;
  }

  template <class F1, class F2, class F3>
  void _initialize_stored_vineyard(int numberOfUpdates, F1&& get_birth, F2&& get_death, F3&& store_cycle) {
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
        vineyard_[bar.dim].push_back({std::forward<F1>(get_birth)(bar.birth), std::forward<F2>(get_death)(bar.death)});
        if (std::forward<F3>(store_cycle)(bar)) {
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
        vineyard_.back().add_pair(std::forward<F1>(get_birth)(bar.birth), std::forward<F2>(get_death)(bar.death));
        if (bar.dim >= static_cast<Dimension>(numberOfBars_.size())) numberOfBars_.resize(bar.dim + 1, 0);
        ++numberOfBars_[bar.dim];
        if (std::forward<F3>(store_cycle)(bar)) {
          auto cycle = base_.get_current_representative_cycle(idx, true);
          latest_representative_cycles_->emplace_back(Cycle(cycle.begin(), cycle.end()), bar.dim, idx);
        }
        ++idx;
      }
    }
  }

  template <class F1, class F2, class F3>
  void _update_stored_vineyard(F1&& get_birth, F2&& get_death, F3&& store_cycle) {
    const auto& barcode = base_.get_current_barcode();  // forward only range + order is preserved
    Index idx = 0;
    std::vector<Index> latestIndex;
    if constexpr (flat) latestIndex.resize(vineyard_.size(), -1);
    if (latest_representative_cycles_) {
      latest_representative_cycles_->clear();
    }
    for (const auto& bar : barcode) {
      T birth = std::forward<F1>(get_birth)(bar.birth);
      T death = std::forward<F2>(get_death)(bar.death);
      if constexpr (flat) {
        ++latestIndex[bar.dim];
        vineyard_[bar.dim].push_back({birth, death});
      } else {
        vineyard_[idx].add_pair(birth, death);
      }
      if (std::forward<F3>(store_cycle)(bar)) {
        const auto& cycle = base_.get_current_representative_cycle(idx, true);
        latest_representative_cycles_->emplace_back(Cycle(cycle.begin(), cycle.end()), bar.dim,
                                                    flat ? latestIndex[bar.dim] : idx);
      }
      ++idx;
    }
  }
};

}  // namespace vineyard
}  // namespace Gudhi

#endif  // GUDHI_VINEYARD_BUILDER_H_
