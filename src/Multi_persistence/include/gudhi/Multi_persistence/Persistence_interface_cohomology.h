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

#include <boost/range/iterator_range_core.hpp>

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
  using PCOH_complex = Multi_parameter_filtered_complex_pcoh_interface<MultiFiltrationValue>;
  using Complex = typename PCOH_complex::Complex;                                /**< Complex type */
  using Dimension = typename PCOH_complex::Dimension;                            /**< Dimension type */
  using Index = typename PCOH_complex::Simplex_key;                              /**< Index type */
  using Map = typename PCOH_complex::Map;                                        /**< Map type */
  using Bar = Gudhi::persistence_matrix::Persistence_interval<Dimension, Index>; /**< Bar type */
  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
  using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<PCOH_complex, Field_Zp>;
  using Barcode = std::vector<Bar>;                                              /**< Barcode type */
  template<class Complex>
  using As_type = Persistence_interface_cohomology<typename Complex::Filtration_value>; /**< This type. */

  static constexpr const auto nullDeath = Bar::inf;
  static constexpr const bool is_vine = false;        /** False. */
  static constexpr const bool has_rep_cycles = false; /** False. */

  Persistence_interface_cohomology() : interface_(), barcode_() {}

  // `permutation` is assumed to have stable size, i.e., its address never changes
  Persistence_interface_cohomology(const Complex& cpx, const Map& permutation) : interface_(cpx, permutation)
  {
    _initialize();
  }

  Persistence_interface_cohomology(const Persistence_interface_cohomology& other) = delete;

  // permutation is assumed to be the same than from the copied object, just its address can change
  Persistence_interface_cohomology(const Persistence_interface_cohomology& other, const Map& permutation)
      : interface_(other.interface_, permutation), barcode_(other.barcode_)
  {}

  Persistence_interface_cohomology(Persistence_interface_cohomology&& other) = delete;

  // permutation is assumed to be the same than from the moved object, just its address can change
  Persistence_interface_cohomology(Persistence_interface_cohomology&& other, const Map& permutation)
      : interface_(std::move(other.interface_), permutation), barcode_(std::move(other.barcode_))
  {}

  ~Persistence_interface_cohomology() = default;

  Persistence_interface_cohomology& operator=(const Persistence_interface_cohomology& other) = delete;
  Persistence_interface_cohomology& operator=(Persistence_interface_cohomology&& other) noexcept = delete;

  // TODO: swap?

  template <class Complex>
  void reinitialize(const Complex& cpx, const Map& permutation)
  {
    interface_.reinitialize(cpx, permutation);
    _initialize();
  }

  void reset()
  {
    interface_.reset();
    barcode_.clear();
  }

  [[nodiscard]] bool is_initialized() const { return interface_.is_initialized(); }

  Dimension get_dimension(Index i) const
  {
    GUDHI_CHECK(is_initialized(), "Dimension can not be computed uninitialized.");
    return interface_.dimension(interface_.simplex(i));
  }

  const Barcode& get_barcode()
  {
    GUDHI_CHECK(is_initialized(), "Barcode can not be computed uninitialized.");
    return barcode_;
  }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream& operator<<(std::ostream& stream, const Persistence_interface_cohomology& pers)
  {
    stream << "Complex:\n";
    stream << pers.interface_ << "\n";
    stream << "Barcode:\n";
    for (const auto bar : pers.barcode_) {
      stream << bar << "\n";
    }
    stream << "\n";

    return stream;
  }

 private:
  PCOH_complex interface_;
  Barcode barcode_;

  void _initialize()
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
