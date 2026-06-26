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
 * @file Persistence_interface_homology.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::multi_persistence::Persistence_interface_homology class.
 */

#ifndef MP_PERSISTENCE_INTERFACE_HOMOLOGY_H_INCLUDED
#define MP_PERSISTENCE_INTERFACE_HOMOLOGY_H_INCLUDED

#include <numeric>
#include <utility>
#include <vector>

#include <boost/range/adaptor/transformed.hpp>

#include <gudhi/Matrix.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Persistence_interface_homology Persistence_interface_homology.h \
 * gudhi/Multi_persistence/Persistence_interface_homology.h
 * @ingroup multi_persistence
 *
 * @brief Interface respecting the @ref PersistenceAlgorithm concept to use @ref Slicer with the homology
 * and representative cycle algorithms implemented in @ref Gudhi::persistence_matrix::Matrix.
 *
 * @tparam PosIdxPersistenceMatrixOptions Options respecting the
 * @ref Gudhi::persistence_matrix::PersistenceMatrixOptions concept such that
 * @ref Gudhi::persistence_matrix::PersistenceMatrixOptions::has_column_pairings is true. Is it also recommended
 * that @ref Gudhi::persistence_matrix::PersistenceMatrixOptions::Index and
 * @ref Gudhi::persistence_matrix::PersistenceMatrixOptions::Dimension correspond to the respective types in the used
 * @ref Slicer.
 * @tparam MultiFiltrationValue Value of @ref Slicer::Filtration_value.
 */
template <class PosIdxPersistenceMatrixOptions, class MultiFiltrationValue>
class Persistence_interface_homology
{
 public:
  using Options = PosIdxPersistenceMatrixOptions;            /**< Matrix options */
  using Dimension = typename Options::Dimension;             /**< Dimension type */
  using Index = typename Options::Index;                     /**< Index type */
  using Complex = Multi_parameter_filtered_complex<MultiFiltrationValue, Index, Dimension>; /**< Complex type */
  using Matrix = Gudhi::persistence_matrix::Matrix<Options>; /**< Matrix type */
  using Map = std::vector<Index>;                            /**< Permutation map for filtration order */
  using Bar = typename Matrix::Bar;                          /**< Bar type */
  using Cycle = typename Matrix::Cycle;                      /**< Cycle type */
  template <class Other_complex>
  using As_type = Persistence_interface_homology<PosIdxPersistenceMatrixOptions,
                                                 typename Other_complex::Filtration_value>; /**< This type. */

  static constexpr auto nullDeath = Bar::inf;               /**< Death value of the barcode when the bar never died. */
  static constexpr auto nullDimension = Matrix::template get_null_value<Dimension>(); /**< None value for dimension. */
  static constexpr bool is_vine = false;                                              /**< False: no cheap update. */
  /**
   * @brief True if and only if PosIdxPersistenceMatrixOptions::can_retrieve_representative_cycles is true.
   */
  static constexpr bool has_rep_cycles = Options::can_retrieve_representative_cycles;

  Persistence_interface_homology() : complex_(nullptr) {}

  template <class Filtration_range>
  Persistence_interface_homology(const Complex& cpx, const Filtration_range& filtrationValues, bool ignoreInf = false)
      : complex_(&cpx), matrix_(filtrationValues.size()), order_(filtrationValues.size())
  {
    static_assert(Options::has_column_pairings, "Underlying matrix has to store barcode.");

    _initialize_order(filtrationValues, ignoreInf);
    _initialize_persistence();
  }

  Persistence_interface_homology(const Persistence_interface_homology& other) = default;

  Persistence_interface_homology(Persistence_interface_homology&& other) noexcept
      : complex_(std::exchange(other.complex_, nullptr)),
        matrix_(std::move(other.matrix_)),
        order_(std::move(other.order_))
  {}

  ~Persistence_interface_homology() = default;

  Persistence_interface_homology& operator=(const Persistence_interface_homology& other) = default;

  Persistence_interface_homology& operator=(Persistence_interface_homology&& other) noexcept
  {
    complex_ = std::exchange(other.complex_, nullptr);
    matrix_ = std::move(other.matrix_);
    order_ = std::move(other.order_);
    return *this;
  }

  // TODO: swap?

  /** 
   * @brief The `ignoreInf` argument is not ignored for this class.
   */
  template <class Filtration_range>
  void initialize(const Complex& cpx, const Filtration_range& filtrationValues, bool ignoreInf = false)
  {
    complex_ = &cpx;
    matrix_ = Matrix(filtrationValues.size());
    order_.resize(filtrationValues.size());

    _initialize_order(filtrationValues, ignoreInf);
    _initialize_persistence();
  }

  /** 
   * @brief The `ignoreInf` argument is not ignored for this class.
   */
  template <class Filtration_range>
  void update(const Filtration_range& filtrationValues, bool ignoreInf = false)
  {
    GUDHI_CHECK(is_initialized(), std::logic_error("Barcode can not be updated uninitialized."));

    initialize(*complex_, filtrationValues, ignoreInf);
  }

  [[nodiscard]] bool is_initialized() const { return complex_ != nullptr; }

  const Map& get_current_order() const { return order_; }

  auto get_barcode()
  {
    GUDHI_CHECK(is_initialized(), std::logic_error("Barcode can not be computed uninitialized."));

    const auto& barcode = matrix_.get_current_barcode();
    return boost::adaptors::transform(barcode, [&](const Bar& bar) -> Bar {
      return Bar(order_[bar.birth], bar.death == Bar::inf ? bar.death : order_[bar.death], bar.dim);
    });
  }

  auto get_all_representative_cycles(bool update = true, Dimension dim = nullDimension)
  {
    static_assert(has_rep_cycles, "`get_all_representative_cycles` is not enabled with the given options.");
    GUDHI_CHECK(is_initialized(), std::logic_error("Representative cycles can not be computed uninitialized."));

    if (update) matrix_.update_all_representative_cycles(dim);
    const auto& cycles = matrix_.get_all_representative_cycles();
    return boost::adaptors::transform(cycles, [&](const Cycle& cycle) -> Cycle {
      Cycle c(cycle.size());
      for (Index i = 0; i < cycle.size(); ++i) {
        c[i] = order_[cycle[i]];
      }
      return c;
    });
  }

  auto get_representative_cycle(Index barcodeIndex, bool update = true)
  {
    static_assert(has_rep_cycles, "`get_representative_cycle` is not enabled with the given options.");
    GUDHI_CHECK(is_initialized(), std::logic_error("Representative cycles can not be computed uninitialized."));

    const auto& bar = matrix_.get_current_barcode()[barcodeIndex];
    if (update) matrix_.update_representative_cycle(bar);
    const auto& cycle = matrix_.get_representative_cycle(bar);
    return boost::adaptors::transform(cycle, [&](const Index& i) -> Index { return order_[i]; });
  }

  friend std::ostream& operator<<(std::ostream& stream, Persistence_interface_homology& pers)
  {
    stream << "Matrix:\n";
    stream << "[\n";
    for (auto i = 0U; i < pers.matrix_.get_number_of_columns(); i++) {
      stream << "[";
      for (const auto& v : pers.matrix_.get_column(i)) stream << v << ", ";
      stream << "]\n";
    }
    stream << "]\n";
    stream << "Permutation:\n";
    for (auto v : *pers.order_) {
      stream << v << " ";
    }
    stream << "\n";

    return stream;
  }

 private:
  Complex const* complex_; /**< Pointer to complex. */
  Matrix matrix_;
  Map order_;

  template <class Filtration_range>
  void _initialize_order(const Filtration_range& filtrationValues, bool ignoreInf)
  {
    const auto& dimensions = complex_->get_dimensions();

    GUDHI_CHECK(complex_->get_boundaries().size() == dimensions.size(),
                std::invalid_argument("Boundary and dimension range sizes are not matching."));
    GUDHI_CHECK(complex_->get_boundaries().size() == filtrationValues.size(),
                std::invalid_argument("Boundary and filtration value range sizes are not matching."));

    std::iota(order_.begin(), order_.end(), 0);
    std::sort(order_.begin(), order_.end(), [&](Index i, Index j) {
      if (ignoreInf) {
        if (filtrationValues[i] != MultiFiltrationValue::T_inf && filtrationValues[j] == MultiFiltrationValue::T_inf)
          return true;
        // all elements at inf are considered equal
        if (filtrationValues[i] == MultiFiltrationValue::T_inf) return false;
      }
      if (dimensions[i] < dimensions[j]) return true;
      if (dimensions[i] > dimensions[j]) return false;
      return filtrationValues[i] < filtrationValues[j];
    });
    if (ignoreInf) {
      Index end = order_.size();
      while (end > 0 && filtrationValues[order_[end - 1]] == MultiFiltrationValue::T_inf) --end;
      order_.resize(end);
    }
  }

  void _initialize_persistence()
  {
    const auto& boundaryMatrix = complex_->get_boundaries();
    const auto& dimensions = complex_->get_dimensions();

    // simplex IDs need to be increasing in order, so the original ones cannot be used
    Map orderInv(boundaryMatrix.size());
    Map translatedBoundary;
    Index id = 0;
    for (auto i : order_) {
      orderInv[i] = id;  // order is assumed to be a valid filtration
      translatedBoundary.resize(boundaryMatrix[i].size());
      for (Index j = 0; j < boundaryMatrix[i].size(); ++j) {
        translatedBoundary[j] = orderInv[boundaryMatrix[i][j]];
      }
      std::sort(translatedBoundary.begin(), translatedBoundary.end());
      matrix_.insert_boundary(id, translatedBoundary, dimensions[i]);
      ++id;  // IDs corresponds to the indices in order_
    }
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_PERSISTENCE_INTERFACE_HOMOLOGY_H_INCLUDED
