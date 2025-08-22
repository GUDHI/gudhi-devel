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

/**
 * @brief Backend class computing persistence for @ref Gudhi::multi_persistence::Slicer and eventually vineyards and
 * representative cycles. The concept is realized by @ref Persistence_interface_cohomology and
 * @ref Persistence_interface_matrix.
 */
class PersistenceAlgorithm
{
 public:
  using Dimension = undefined; /**< Dimension type. */
  using Index = undefined;     /**< Index type. */
  using Map = undefined;       /**< Map type. */
  using Bar = undefined;       /**< Bar type. Has to have public members `dim`, `birth` and `death`. */
  template<class Complex>
  using As_type = undefined;   /**< Type of the class for the given template Complex type to be used. */

  /**
   * @brief Barcode type. Has to have a `begin`, `end` and `size` method. When dereferenced, the iterators have
   * to return a object of type @ref Bar.
   */
  using Barcode = undefined;
  /**
   * @brief Cycle type. Only necessary if  @ref has_rep_cycles is true. Has to have an `empty` and `operator[]` method
   * and a copy constructor. The `operator[](i)` has to return the index in the originally given complex of the
   * \f$ i^{th} \f$ element in the cycle.
   */
  using Cycle = undefined;
  /**
   * @brief Cycle container type. Only necessary if  @ref has_rep_cycles is true. Has to have a `begin`, `end` and
   * `size` method. When dereferenced, the iterators have to return a object of type @ref Cycle.
   */
  using Cycles = undefined;

  static constexpr const auto nullDeath = undefined;      /**< Death value of the barcode when the bar never died. */
  static constexpr const bool is_vine = undefined;        /**< True if and only if @ref vine_swap is enabled. */
  static constexpr const bool has_rep_cycles = undefined; /**< True if and only if @ref get_representative_cycles
                                                               is enabled. */

  /**
   * @brief Default constructor. Should be as empty as possible. Except @ref is_initialized (which returns false here),
   * all other method can have an undefined behaviour.
   */
  PersistenceAlgorithm();

  /**
   * @brief Real constructor. Will store a pointer to the permutation map and eventually to the complex.
   * So the addresses are assumed to remain valid until the end of this object. If they change, they can be updated
   * via the copy/move constructors if the content did not change and through @ref reinitialize otherwise.
   * 
   * @tparam Complex Any version of @ref Multi_parameter_filtered_complex.
   * @param cpx Complex containing boundaries and dimensions.
   * @param permutation Permutation map indicating the order of the cells stored in @p cpx as a standard
   * 1-dimensional filtration. I.e., `permutation[i]` corresponds to the index in `cpx` which should be
   * used as the i-th cell in the filtration.
   */
  template <class Complex>
  PersistenceAlgorithm(const Complex& cpx, const Map& permutation);

  /**
   * @brief Standard copy constructor deleted to avoid accidents with invalidated addresses.
   */
  PersistenceAlgorithm(const PersistenceAlgorithm& other) = delete;

  /**
   * @brief Copy constructor. Can assumes that the complex used at construction did not change address and updates the
   * permutation pointer with the given map if the complex was initialized.
   *
   * @param other To copy.
   * @param permutation Permutation map. Has to correspond to the same map than in @p other, except that its address
   * does not have to be the same.
   */
  PersistenceAlgorithm(const PersistenceAlgorithm& other, const Map& permutation);

  /**
   * @brief Standard move constructor deleted to avoid accidents with invalidated addresses.
   */
  PersistenceAlgorithm(PersistenceAlgorithm&& other) = delete;

  /**
   * @brief Move constructor. Can assumes that the complex used at construction did not change address and updates the
   * permutation pointer with the given map if the complex was initialized.
   *
   * @param other To move.
   * @param permutation Permutation map. Has to correspond to the same map than in @p other, except that its address
   * does not have to be the same.
   */
  PersistenceAlgorithm(PersistenceAlgorithm&& other, const Map& permutation);

  // TODO: swap?

  /**
   * @brief Reinitializes the interface with the new given complex and permutation.
   * To use instead of the classical assign operator `operator=`.
   */
  template <class Complex>
  void reinitialize(const Complex& cpx, const Map& permutation);

  /**
   * @brief Empties the class such that @ref is_initialized returns false.
   */
  void reset();

  /**
   * @brief Returns `true` if and only if all stored pointers are not null.
   */
  bool is_initialized() const;

  /**
   * @brief Returns the dimension of the \f$ i^{th} \f$ cell in the filtration.
   */
  Dimension get_dimension(Index i) const;

  /**
   * @brief Returns the barcode of the filtration.
   */
  Barcode get_barcode();

  /**
   * @brief Only enabled if @ref is_vine is true. Swaps the \f$ i^{th} \f$ and \f$ i^{th} + 1 \f$ cell's position in
   * the filtration while updating the barcode.
   */
  void vine_swap(Index i);

  /**
   * @brief Only enabled if @ref has_rep_cycles is true. Returns the representative cycles of the current barcode.
   * 
   * @param update If true, updates the stored representative cycles, otherwise just returns the container in its
   * current state. So should be true at least the first time the method is used.
   */
  Cycles get_representative_cycles(bool update);

  /**
   * @brief Only necessary if the `operator<<` method of the Slicer is needed.
   */
  friend std::ostream& operator<<(std::ostream& stream, PersistenceAlgorithm& pers);
};

}  // namespace multi_persistence
}  // namespace Gudhi
