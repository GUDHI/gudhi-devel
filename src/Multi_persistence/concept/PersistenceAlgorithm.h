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
 * representative cycles. The concept is realized by @ref Persistence_interface_cohomology,
 * @ref Persistence_interface_homology and @ref Persistence_interface_vineyard.
 */
class PersistenceAlgorithm
{
 public:
  using Dimension = undefined; /**< Dimension type. */
  using Index = undefined;     /**< Index type. */
  using Map = undefined;       /**< Map type. */
  template <class Complex>
  using As_type = undefined;   /**< Type of this class with given template Complex type. */

  static constexpr auto nullDeath = undefined;      /**< Death value of the barcode when the bar never died. */
  static constexpr bool has_rep_cycles = undefined; /**< True if and only if rep. cycle methods are enabled. */

  /**
   * @brief Default constructor. Should be as empty as possible. Except @ref is_initialized (which returns false here),
   * all other method can have an undefined behaviour.
   */
  PersistenceAlgorithm();

  /**
   * @brief Standard copy constructor.
   */
  PersistenceAlgorithm(const PersistenceAlgorithm& other);

  /**
   * @brief Standard move constructor.
   */
  PersistenceAlgorithm(PersistenceAlgorithm&& other) noexcept;

  /**
   * @brief Standard destructor
   */
  ~PersistenceAlgorithm();

  /**
   * @brief Standard copy assign
   */
  PersistenceAlgorithm& operator=(const PersistenceAlgorithm& other);

  /**
   * @brief Standard move assign
   */
  PersistenceAlgorithm& operator=(PersistenceAlgorithm&& other) noexcept;

  // TODO: swap?

  /**
   * @brief Initializes the persistence computation with the given complex and filtration values. The filtration values
   * have to represent a valid 1-dimensional filtration, otherwise the behaviour is undefined.
   *
   * @tparam Complex @ref Multi_parameter_filtered_complex with desired template parameters
   * @tparam Filtration_range Range of Complex::Filtration_value::value_type with `size` and `operator[]` method.
   * @param cpx Complex containing all cells of the filtration.
   * @param filtrationValues Range of filtration values with indices corresponding to the indices in @p cpx.
   * @param ignoreInf If true, all cells at infinity filtration values are ignored for the initialization. Can be
   * unused in the method if ignoring values is not appropriate for the update method.
   */
  template <class Complex, class Filtration_range>
  void initialize(const Complex& cpx, const Filtration_range& filtrationValues, [[maybe_unused]] bool ignoreInf);

  /**
   * @brief Updates the persistence barcode with the given filtration values, assuming the complex remains the same
   * than for the initialization.
   *
   * @tparam Filtration_range Range of Complex::Filtration_value::value_type (used for @ref initialize) with `size`
   * and `operator[]` method.
   * @param filtrationValues Range of filtration values with indices corresponding to the indices in originally given
   * complex.
   * @param ignoreInf If true, all cells at infinity filtration values are ignored for the initialization. Can be
   * unused in the method if ignoring values is not appropriate for the update method.
   */
  template <class Filtration_range>
  void update(const Filtration_range& filtrationValues, [[maybe_unused]] bool ignoreInf);

  /**
   * @brief Returns `true` if and only if @ref initialize was called at least once.
   */
  bool is_initialized() const;

  /**
   * @brief Returns the current permutation map indicating the filtration order. Has a `size` and `operator[]` method.
   */
  const Map& get_current_order() const;

  /**
   * @brief Returns the barcode of the filtration. The barcode must have a `begin`, `end` and `size` method. A bar in
   * the barcode has to be of type a variation of @ref Gudhi::persistence_matrix::Persistence_interval.
   */
  auto get_barcode();

  /**
   * @brief Only enabled if @ref has_rep_cycles is true. Returns the representative cycles of the current barcode.
   *
   * @param update If true, updates the stored representative cycles, otherwise just returns the container in its
   * current state. So should be true at least the first time the method is used.
   * @param dim If given, returns only the cycles in this dimension. Otherwise, returns all cycles.
   */
  auto get_all_representative_cycles(bool update, Dimension dim = nullDimension);

  /**
   * @brief Only enabled if @ref has_rep_cycles is true. Returns the representative cycle at given index in the current
   * barcode.
   *
   * @param barcodeIndex Index of the bar in @ref get_barcode.
   * @param update If true, updates the stored representative cycles, otherwise just returns the container in its
   * current state. So should be true at least the first time the method is used.
   */
  auto get_representative_cycle(Index barcodeIndex, bool update);

  /**
   * @brief Only necessary if the `operator<<` method of the Slicer is needed.
   */
  friend std::ostream& operator<<(std::ostream& stream, PersistenceAlgorithm& pers);
};

}  // namespace multi_persistence
}  // namespace Gudhi
