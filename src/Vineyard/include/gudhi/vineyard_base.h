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
 * @file vineyard_base.h
 * @author Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::vineyard::Vineyard_base class.
 */

#ifndef GUDHI_VINEYARD_BASE_H_
#define GUDHI_VINEYARD_BASE_H_

#include <numeric>
#include <vector>

#include <boost/range/adaptor/transformed.hpp>

#include <gudhi/Debug_utils.h>
#include <gudhi/Matrix.h>

namespace Gudhi {
namespace vineyard {

/**
 * @ingroup vineyard
 *
 * @brief Options for the internal matrix of @ref Vineyard_base.
 *
 * @tparam column_type Column type of the matrix.
 */
template <typename D, typename I, bool is_RU, Gudhi::persistence_matrix::Column_types column_type>
struct Vineyard_matrix_options : Gudhi::persistence_matrix::Default_options<column_type, true> {
  using Dimension = D;
  using Index = I;

  static const bool is_of_boundary_type = is_RU;
  static const Gudhi::persistence_matrix::Column_indexation_types column_indexation_type =
      Gudhi::persistence_matrix::Column_indexation_types::POSITION;
  static const bool has_column_pairings = true;
  static const bool has_vine_update = true;
  static const bool can_retrieve_representative_cycles = true;
};

/**
 * @ingroup vineyard
 *
 * @brief Default options for @ref Vineyard_base.
 */
struct Default_vineyard_options {
  using Index = std::uint32_t; /**< Index type. */
  using Dimension = int;       /**< Dimension value type. */

  static constexpr bool is_RU = false;  // TODO: benchmark both
  /**
   * @brief Column type use by the internal matrix.
   */
  static const Gudhi::persistence_matrix::Column_types column_type =
      Gudhi::persistence_matrix::Column_types::NAIVE_VECTOR;  // TODO: benchmark
};

/**
 * @class Vineyard_base vineyard_base.h gudhi/vineyard_base.h
 * @brief
 * @details
 *
 * @ingroup vineyard
 *
 * @tparam VineyardOptions Structure following the @ref VineyardOptions concept. Default value: @ref
 * Default_vineyard_options.
 */
template <class VineyardOptions = Default_vineyard_options>
class Vineyard_base
{
 public:
  using Index = typename VineyardOptions::Index;         /**< Complex index type. */
  using Dimension = typename VineyardOptions::Dimension; /**< Dimension type. */
  using Matrix_options =
      Vineyard_matrix_options<Dimension, Index, VineyardOptions::is_RU, VineyardOptions::column_type>;
  using Matrix = Gudhi::persistence_matrix::Matrix<Matrix_options>;
  using Bar = typename Matrix::Bar;     /**< Bar type. */
  using Cycle = typename Matrix::Cycle; /**< Cycle type. */
  using Permutation = std::vector<Index>;

  static constexpr Dimension nullDimension = Matrix::template get_null_value<Dimension>();

  Vineyard_base() = default;

  template <class Boundary_range, class Dimension_range, class Filtration_range>
  Vineyard_base(const Boundary_range& boundaryMatrix,
                const Dimension_range& dimensions,
                const Filtration_range& filtrationValues)
      : matrix_(boundaryMatrix.size()), order_(boundaryMatrix.size())
  {
    // All static_assert in this class are quite useless as Matrix_options is fixed and has those enabled
    // I keep them just here for now in case Matrix_options becomes a template argument instead later
    static_assert(Matrix_options::has_vine_update, "Underlying matrix has to support vine swaps.");
    static_assert(Matrix_options::has_column_pairings, "Underlying matrix has to store barcode.");
    static_assert(
        Matrix_options::column_indexation_type == Gudhi::persistence_matrix::Column_indexation_types::POSITION ||
            (Matrix_options::is_of_boundary_type &&
             Matrix_options::column_indexation_type == Gudhi::persistence_matrix::Column_indexation_types::CONTAINER),
        "Matrix has a non supported index scheme.");

    if constexpr (!Matrix_options::is_of_boundary_type) {
      idToPos_.emplace();
      idToPos_->resize(order_.size());
    }

    _initialize(boundaryMatrix, dimensions, filtrationValues);
  }

  template <class Boundary_range, class Dimension_range, class Filtration_range>
  void initialize(const Boundary_range& boundaryMatrix,
                  const Dimension_range& dimensions,
                  const Filtration_range& filtrationValues)
  {
    matrix_ = Matrix(boundaryMatrix.size());
    order_.resize(boundaryMatrix.size());
    if constexpr (!Matrix_options::is_of_boundary_type) {
      if (!idToPos_) idToPos_.emplace();
      idToPos_->resize(order_.size());
    }

    _initialize(boundaryMatrix, dimensions, filtrationValues);
  }

  template <class Filtration_range>
  void update(const Filtration_range& filtrationValues)
  {
    GUDHI_CHECK(filtrationValues.size() == order_.size(),
                std::invalid_argument("Filtration value container size is not matching."));

    for (Index i = 1; i < order_.size(); i++) {
      int curr = i;
      // speed up when ordered by dim, to avoid unnecessary swaps
      while (curr > 0 && matrix_.get_column_dimension(curr) == matrix_.get_column_dimension(curr - 1) &&
             filtrationValues[order_[curr]] < filtrationValues[order_[curr - 1]]) {
        if constexpr (!Matrix_options::is_of_boundary_type) {
          auto id1 = matrix_.get_pivot(curr - 1);
          auto id2 = matrix_.get_pivot(curr);
          std::swap((*idToPos_)[id1], (*idToPos_)[id2]);
        }
        matrix_.vine_swap(curr - 1);
        std::swap(order_[curr - 1], order_[curr]);
        --curr;
      }
    }
  }

  [[nodiscard]] bool is_initialized() const { return !order_.empty(); }

  const Permutation& get_current_order() const { return order_; }

  auto get_current_barcode()
  {
    const auto& barcode = matrix_.get_current_barcode();
    return boost::adaptors::transform(barcode, [&](const Bar& bar) -> Bar {
      return Bar(order_[bar.birth], bar.death == Bar::inf ? bar.death : order_[bar.death], bar.dim);
    });
  }

  auto get_all_current_representative_cycles(bool update = true, Dimension dim = nullDimension)
  {
    static_assert(Matrix_options::can_retrieve_representative_cycles,
                  "Underlying matrix has to support representative cycles.");

    if (update) matrix_.update_all_representative_cycles(dim);
    const auto& cycles = matrix_.get_all_representative_cycles();
    return boost::adaptors::transform(cycles, [&](const Cycle& cycle) -> Cycle {
      Cycle c(cycle.size());
      for (Index i = 0; i < cycle.size(); ++i) {
        if constexpr (Matrix_options::is_of_boundary_type) {
          // works for RU because id == pos, but does not work for chain with vine
          // we need a id to pos map in that case
          c[i] = order_[cycle[i]];
        } else {
          c[i] = order_[(*idToPos_)[cycle[i]]];
        }
      }
      return c;
    });
  }

  auto get_current_representative_cycle(Index barcodeIndex, bool update = true)
  {
    static_assert(Matrix_options::can_retrieve_representative_cycles,
                  "Underlying matrix has to support representative cycles.");

    const auto& bar = matrix_.get_current_barcode()[barcodeIndex];
    if (update) matrix_.update_representative_cycle(bar);
    const auto& cycle = matrix_.get_representative_cycle(bar);
    return boost::adaptors::transform(cycle, [&](const Index& i) -> Index {
      if constexpr (Matrix_options::is_of_boundary_type) {
        // works for RU because id == pos, but does not work for chain with vine
        // we need a id to pos map in that case
        return order_[i];
      } else {
        return order_[(*idToPos_)[i]];
      }
    });
  }

  // debug purposes mainly
  /**
   * @brief Outstream operator.
   */
  friend std::ostream& operator<<(std::ostream& stream, Vineyard_base& vyd)
  {
    stream << "Matrix:\n";
    stream << "[\n";
    for (auto i = 0U; i < vyd.matrix_.get_number_of_columns(); i++) {
      stream << "[";
      for (const auto& v : vyd.matrix_.get_column(i)) stream << v << ", ";
      stream << "]\n";
    }
    stream << "]\n";
    stream << "Permutation:\n";
    for (auto v : vyd.order_) {
      stream << v << " ";
    }
    stream << "\n";
    if constexpr (Matrix_options::is_of_boundary_type) {
      stream << "ID to position map:\n";
      if (vyd.idToPos_) {
        for (auto v : *vyd.idToPos_) {
          stream << v << " ";
        }
      }
      stream << "\n";
    }
    return stream;
  }

 private:
  Matrix matrix_;
  Permutation order_;
  std::optional<Permutation> idToPos_;  // TODO: remove if chain does not improve run times in benchmark

  template <class Boundary_range, class Dimension_range, class Filtration_range>
  void _initialize(const Boundary_range& boundaryMatrix,
                   const Dimension_range& dimensions,
                   const Filtration_range& filtrationValues)
  {
    GUDHI_CHECK(boundaryMatrix.size() == dimensions.size(),
                std::invalid_argument("Boundary and dimension range sizes are not matching."));
    GUDHI_CHECK(boundaryMatrix.size() == filtrationValues.size(),
                std::invalid_argument("Boundary and filtration value range sizes are not matching."));

    std::iota(order_.begin(), order_.end(), 0);
    std::sort(order_.begin(), order_.end(), [&](Index i, Index j) {
      if (dimensions[i] < dimensions[j]) return true;
      if (dimensions[i] > dimensions[j]) return false;
      return filtrationValues[i] < filtrationValues[j];
    });

    // simplex IDs need to be increasing in order, so the original ones cannot be used
    Permutation orderInv(boundaryMatrix.size());
    Permutation translatedBoundary;
    Index id = 0;
    for (auto i : order_) {
      orderInv[i] = id;  // order is assumed to be a valid filtration
      translatedBoundary.resize(boundaryMatrix[i].size());
      for (Index j = 0; j < boundaryMatrix[i].size(); ++j) {
        translatedBoundary[j] = orderInv[boundaryMatrix[i][j]];
      }
      std::sort(translatedBoundary.begin(), translatedBoundary.end());
      matrix_.insert_boundary(id, translatedBoundary, dimensions[i]);
      if constexpr (!Matrix_options::is_of_boundary_type) {
        (*idToPos_)[id] = id; // before any vine swaps, id == pos
      }
      ++id;  // IDs corresponds to the indices in order_
    }
  }
};

}  // namespace vineyard
}  // namespace Gudhi

#endif  // GUDHI_VINEYARD_BASE_H_
