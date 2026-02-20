/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Persistence_interface_vineyard.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::multi_persistence::Persistence_interface_vineyard class.
 */

#ifndef MP_PERSISTENCE_INTERFACE_VINEYARD_H_INCLUDED
#define MP_PERSISTENCE_INTERFACE_VINEYARD_H_INCLUDED

#include <stdexcept>

#include <gudhi/Debug_utils.h>
#include <gudhi/vineyard_base.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Persistence_interface_vineyard Persistence_interface_vineyard.h \
 * gudhi/Multi_persistence/Persistence_interface_vineyard.h
 * @ingroup multi_persistence
 *
 * @brief Interface respecting the @ref PersistenceAlgorithm concept to use @ref Slicer with the homology, vineyard
 * and representative cycle algorithms implemented in @ref Gudhi::persistence_matrix::Matrix.
 *
 * @tparam VineyardOptions Options respecting the @ref VineyardOptions concept.
 */
template <class VineyardOptions = Gudhi::vineyard::Default_vineyard_options>
class Persistence_interface_vineyard
{
 public:
  using Options = VineyardOptions;
  using Vineyard = Gudhi::vineyard::Vineyard_base<Options>; /**< Vineyard computation base */
  using Dimension = typename Vineyard::Dimension;           /**< Dimension type */
  using Index = typename Vineyard::Index;                   /**< Index type */
  using Bar = typename Vineyard::Bar;                       /**< Bar type */
  using Cycle = typename Vineyard::Cycle;                   /**< Cycle type */
  using Map = typename Vineyard::Permutation;               /**< Permutation map for filtration order. */
  template <class Complex>
  using As_type = Persistence_interface_vineyard<Options>;  /**< This type. */

  static constexpr auto nullDeath = Bar::inf;  /**< Death value of the barcode when the bar never died. */
  static constexpr bool is_vine = true;        /**< True: update optimization + barcode matching enabled. */
  static constexpr bool has_rep_cycles = true; /**< True: rep. cycle methods enabled. */

  Persistence_interface_vineyard() = default;

  template <class Complex, class Filtration_range>
  Persistence_interface_vineyard(const Complex& cpx,
                                 const Filtration_range& filtrationValues,
                                 [[maybe_unused]] bool ignoreInf = false)
      : vineyard_(cpx.get_boundaries(), cpx.get_dimensions(), filtrationValues)
  {}

  Persistence_interface_vineyard(const Persistence_interface_vineyard& other) = default;

  Persistence_interface_vineyard(Persistence_interface_vineyard&& other) noexcept = default;

  ~Persistence_interface_vineyard() = default;

  Persistence_interface_vineyard& operator=(const Persistence_interface_vineyard& other) = default;

  Persistence_interface_vineyard& operator=(Persistence_interface_vineyard&& other) noexcept = default;

  // TODO: swap?

  /** 
   * @brief The `ignoreInf` argument is ignored for this class, as it disrupts the matching property of the barcodes
   * after update.
   */
  template <class Complex, class Filtration_range>
  void initialize(const Complex& cpx, const Filtration_range& filtrationValues, [[maybe_unused]] bool ignoreInf = false)
  {
    vineyard_.initialize(cpx.get_boundaries(), cpx.get_dimensions(), filtrationValues);
  }

  /** 
   * @brief The `ignoreInf` argument is ignored for this class, as it disrupts the matching property of the barcodes
   * after update.
   */
  template <class Filtration_range>
  void update(const Filtration_range& filtrationValues, [[maybe_unused]] bool ignoreInf = false)
  {
    GUDHI_CHECK(is_initialized(), std::logic_error("Barcode can not be updated uninitialized."));

    vineyard_.update(filtrationValues);
  }

  [[nodiscard]] bool is_initialized() const { return vineyard_.is_initialized(); }

  const Map& get_current_order() const { return vineyard_.get_current_order(); }

  auto get_barcode()
  {
    GUDHI_CHECK(is_initialized(), std::logic_error("Barcode can not be computed uninitialized."));

    return vineyard_.get_current_barcode();
  }

  auto get_all_representative_cycles(bool update = true, Dimension dim = Vineyard::nullDimension)
  {
    GUDHI_CHECK(is_initialized(), std::logic_error("Representative cycles can not be computed uninitialized."));

    return vineyard_.get_all_current_representative_cycles(update, dim);
  }

  auto get_representative_cycle(Index barcodeIndex, bool update = true)
  {
    GUDHI_CHECK(is_initialized(), std::logic_error("Representative cycles can not be computed uninitialized."));

    return vineyard_.get_current_representative_cycle(barcodeIndex, update);
  }

  friend std::ostream& operator<<(std::ostream& stream, Persistence_interface_vineyard& pers)
  {
    stream << pers.vineyard_;

    return stream;
  }

 private:
  Vineyard vineyard_;
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_PERSISTENCE_INTERFACE_VINEYARD_H_INCLUDED
