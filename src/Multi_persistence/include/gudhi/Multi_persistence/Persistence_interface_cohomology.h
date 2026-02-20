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
 * @file Persistence_interface_cohomology.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::multi_persistence::Persistence_interface_cohomology class.
 */

#ifndef MP_PERSISTENCE_INTERFACE_COHOMOLOGY_H_INCLUDED
#define MP_PERSISTENCE_INTERFACE_COHOMOLOGY_H_INCLUDED

#include <vector>

#include <gudhi/Multi_persistence/Multi_parameter_filtered_complex_pcoh_interface.h>
#include <gudhi/persistence_interval.h>
#include <gudhi/Persistent_cohomology.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Persistence_interface_cohomology Persistence_interface_cohomology.h \
 * gudhi/Multi_persistence/Persistence_interface_cohomology.h
 * @ingroup multi_persistence
 *
 * @brief Interface respecting the @ref PersistenceAlgorithm concept to use @ref Slicer with the cohomology
 * algorithm implemented in @ref Gudhi::persistent_cohomology::Persistent_cohomology.
 *
 * @tparam MultiFiltrationValue Filtration value type used in @ref Slicer.
 */
template <class MultiFiltrationValue>
class Persistence_interface_cohomology
{
 public:
  /**
   * @brief Complex interface for @ref Gudhi::persistent_cohomology::Persistent_cohomology
   */
  using PCOH_complex = Multi_parameter_filtered_complex_pcoh_interface<MultiFiltrationValue>;
  using Complex = typename PCOH_complex::Complex;                                /**< Complex type */
  using Dimension = typename PCOH_complex::Dimension;                            /**< Dimension type */
  using Index = typename PCOH_complex::Simplex_key;                              /**< Index type */
  using Bar = Gudhi::persistence_matrix::Persistence_interval<Dimension, Index>; /**< Bar type */
  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;                       /**< Cohomology field. */
  using Barcode = std::vector<Bar>;                                              /**< Barcode type */
  using Map = std::vector<Index>;                            /**< Permutation map for filtration order. */
  /**
   * @brief Cohomology computation type
   */
  using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<PCOH_complex, Field_Zp>;
  template <class Complex>
  using As_type = Persistence_interface_cohomology<typename Complex::Filtration_value>; /**< This type. */

  static constexpr const auto nullDeath = Bar::inf;   /**< Death value of the barcode when the bar never died. */
  static constexpr const bool is_vine = false;        /**< False: no cheap update. */
  static constexpr const bool has_rep_cycles = false; /**< False: no rep cycle methods enabled. */

  Persistence_interface_cohomology() : order_(), interface_(), barcode_() {}

  template <class Filtration_range>
  Persistence_interface_cohomology(const Complex& cpx, const Filtration_range& filtrationValues, bool ignoreInf = false)
      : order_(cpx.get_number_of_cycle_generators()), interface_(cpx, order_), barcode_()
  {
    _initialize_order(cpx, filtrationValues, ignoreInf);
    _initialize_persistence();
  }

  Persistence_interface_cohomology(const Persistence_interface_cohomology& other)
      : order_(other.order_), interface_(other.interface_, order_), barcode_(other.barcode_)
  {}

  Persistence_interface_cohomology(Persistence_interface_cohomology&& other) noexcept
      : order_(std::move(other.order_)),
        interface_(std::move(other.interface_), order_),
        barcode_(std::move(other.barcode_))
  {}

  ~Persistence_interface_cohomology() = default;

  Persistence_interface_cohomology& operator=(const Persistence_interface_cohomology& other)
  {
    order_ = other.order_;
    interface_ = PCOH_complex(other.interface_, order_);
    barcode_ = other.barcode_;
    return *this;
  }

  Persistence_interface_cohomology& operator=(Persistence_interface_cohomology&& other) noexcept
  {
    order_ = std::move(other.order_);
    interface_ = PCOH_complex(std::move(other.interface_), order_);
    barcode_ = std::move(other.barcode_);
    return *this;
  }

  // TODO: swap?

  /** 
   * @brief The `ignoreInf` argument is not ignored for this class.
   */
  template <class Filtration_range>
  void initialize(const Complex& cpx, const Filtration_range& filtrationValues, bool ignoreInf = false)
  {
    _initialize_order(cpx, filtrationValues, ignoreInf);  // may trigger a relocation for order_
    interface_ = PCOH_complex(cpx, order_);
    _initialize_persistence();
  }

  /**
   * @brief The `ignoreInf` argument is not ignored for this class.
   */
  template <class Filtration_range>
  void update(const Filtration_range& filtrationValues, bool ignoreInf = false)
  {
    GUDHI_CHECK(is_initialized(), std::logic_error("Barcode can not be updated uninitialized."));

    auto const* cpx = interface_.get_complex_ptr();

    // should not trigger a relocation for order_ as complex size did not change
    _initialize_order(*cpx, filtrationValues, ignoreInf);
    _initialize_persistence();
  }

  [[nodiscard]] bool is_initialized() const { return interface_.is_initialized(); }

  const Map& get_current_order() const { return order_; }

  const Barcode& get_barcode()
  {
    GUDHI_CHECK(is_initialized(), std::logic_error("Barcode can not be computed uninitialized."));
    return barcode_;
  }

  friend std::ostream& operator<<(std::ostream& stream, const Persistence_interface_cohomology& pers)
  {
    stream << "Complex:\n";
    stream << pers.interface_ << "\n";
    stream << "Permutation:\n";
    for (auto v : pers.order_) {
      stream << v << " ";
    }
    stream << "\n";
    stream << "Barcode:\n";
    for (const auto bar : pers.barcode_) {
      stream << bar << "\n";
    }
    stream << "\n";

    return stream;
  }

 private:
  Map order_;
  PCOH_complex interface_;
  Barcode barcode_;

  template <class Filtration_range>
  void _initialize_order(const Complex& complex, const Filtration_range& filtrationValues, bool ignoreInf)
  {
    const auto& dimensions = complex.get_dimensions();
    order_.resize(complex.get_number_of_cycle_generators());
    std::iota(order_.begin(), order_.end(), 0);
    std::sort(order_.begin(), order_.end(), [&](Index i, Index j) {
      if (ignoreInf) {
        if (filtrationValues[i] != MultiFiltrationValue::T_inf && filtrationValues[j] == MultiFiltrationValue::T_inf)
          return true;
        // all elements at inf are considered equal
        if (filtrationValues[i] == MultiFiltrationValue::T_inf) return false;
      }
      if (dimensions[i] > dimensions[j]) return false;
      if (dimensions[i] < dimensions[j]) return true;
      // if filtration values are equal, we don't care about order, so considered the same object
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
    Persistent_cohomology pcoh(interface_, true);
    pcoh.init_coefficients(2);
    pcoh.compute_persistent_cohomology_without_optimizations(0);
    const auto& pairs = pcoh.get_persistent_pairs();

    barcode_ = Barcode(pairs.size());
    Index i = 0;
    for (const auto& p : pairs) {
      auto& b = barcode_[i];
      b.dim = interface_.dimension(get<0>(p));
      b.birth = get<0>(p);
      b.death = get<1>(p) == PCOH_complex::null_simplex() ? Bar::inf : get<1>(p);
      ++i;
    }
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_PERSISTENCE_INTERFACE_COHOMOLOGY_H_INCLUDED
