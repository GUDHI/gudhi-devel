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
 * @file PersistenceAlgorithm.h
 * @author Hannah Schreiber
 * @brief Contains the concept for the persistence interface needed for @ref Gudhi::multi_persistence::Slicer.
 */

namespace Gudhi {
namespace multi_persistence {

class PersistenceAlgorithm
{
 public:
  using Dimension = undefined;
  using Index = undefined;
  using Map = undefined;
  using Bar = undefined;
  using Barcode = undefined;
  using Cycle = undefined;
  using Cycles = undefined;

  static constexpr const auto nullDeath = undefined;
  static constexpr const bool is_vine = undefined;
  static constexpr const bool has_rep_cycles = undefined;

  PersistenceAlgorithm();

  // `permutation` is assumed to have stable size, i.e., its address never changes
  template <class Complex>
  PersistenceAlgorithm(const Complex& cpx, const Map& permutation);

  PersistenceAlgorithm(const PersistenceAlgorithm& other) = delete;

  // for thread safe slicer, where cpx is guaranteed to remain its address
  // permutation is assumed to be the same than from the copied object, just its address can change
  PersistenceAlgorithm(const PersistenceAlgorithm& other, const Map& permutation);

  PersistenceAlgorithm(PersistenceAlgorithm&& other) = delete;

  // for thread safe slicer, where cpx is guaranteed to remain its address
  // permutation is assumed to be the same than from the moved object, just its address can change
  PersistenceAlgorithm(PersistenceAlgorithm&& other, const Map& permutation);

  // TODO: swap?

  template <class Complex>
  void reinitialize(const Complex& cpx, const Map& permutation);
  bool is_initialized() const;
  Dimension get_dimension(Index i) const;
  Barcode get_barcode();
  void vine_swap(Index i);
  Cycles get_representative_cycles(bool update);
};

}  // namespace multi_persistence
}  // namespace Gudhi
