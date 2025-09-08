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
 * @file Multi_parameter_filtered_complex_pcoh_interface.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::multi_persistence::Multi_parameter_filtered_complex_pcoh_interface class.
 */

#ifndef MP_COMPLEX_PCOH_INTERFACE_H_INCLUDED
#define MP_COMPLEX_PCOH_INTERFACE_H_INCLUDED

#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_parameter_filtered_complex.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Multi_parameter_filtered_complex_pcoh_interface Multi_parameter_filtered_complex_pcoh_interface.h \
 * gudhi/Multi_persistence/Multi_parameter_filtered_complex_pcoh_interface.h
 * @ingroup multi_persistence
 *
 * @brief Interface respecting the @ref FilteredComplex concept to use
 * @ref Gudhi::persistent_cohomology::Persistent_cohomology with an instantiation of
 * @ref Multi_parameter_filtered_complex.
 *
 * @tparam MultiFiltrationValue Filtration value type used in @ref Multi_parameter_filtered_complex.
 */
template <class MultiFiltrationValue>
class Multi_parameter_filtered_complex_pcoh_interface
{
 public:
  using Complex = Multi_parameter_filtered_complex<MultiFiltrationValue>; /**< Complex type */
  using Simplex_key = std::uint32_t;                                      /**< Simplex_key type */
  using Simplex_handle = Simplex_key;                                     /**< Simplex_handle type */
  using Filtration_value = Simplex_key;                                   /**< Internal filtration value type */
  using Dimension = int;                                                  /**< Internal dimension type */
  using Map = std::vector<Simplex_handle>;                                /**< Map type */
  using Filtration_simplex_range = Map;                                   /**< Filtration_simplex_range type */
  using Boundary_simplex_range = Map;                                     /**< Boundary_simplex_range type */

  /**
   * @brief Default constructor, storing null pointers.Can be tested with @ref is_initialized(), but any other method
   * should not be used with this constructor.
   */
  Multi_parameter_filtered_complex_pcoh_interface() : boundaries_(nullptr), newToOldPerm_(nullptr) {}

  /**
   * @brief Constructs the interface and stores a pointer to the given complex and the given permutation. Assumes that
   * the pointers remain valid through the whole life of the instance. If they get invalidated, they can be updated
   * via the copy/move constructors if the content did not change (otherwise, reconstruct from scratch with
   * @ref reinitialize).
   *
   * @param boundaries Complex to interface.
   * @param permutation Permutation map indicating the order of the cells stored in @p boundaries as a standard
   * 1-dimensional filtration. I.e., `permutation[i]` corresponds to the index in `boundaries` which should be
   * used as the i-th cell in the filtration.
   */
  Multi_parameter_filtered_complex_pcoh_interface(const Complex &boundaries, const Map &permutation)
      : boundaries_(&boundaries), newToOldPerm_(&permutation), keys_(boundaries.get_number_of_cycle_generators(), -1)
  {}

  Multi_parameter_filtered_complex_pcoh_interface(const Multi_parameter_filtered_complex_pcoh_interface &other) =
      delete;

  /**
   * @brief Copy constructor. Copies the complex pointer from @p other and updates the permutation pointer with the
   * given map if the complex was initialized.
   *
   * @param other Interface to copy.
   * @param permutation Permutation map. Has to correspond to the same map than in @p other, except that its address
   * does not have to be the same.
   */
  Multi_parameter_filtered_complex_pcoh_interface(const Multi_parameter_filtered_complex_pcoh_interface &other,
                                                  const Map &permutation)
      : boundaries_(other.boundaries_),
        newToOldPerm_(boundaries_ != nullptr ? &permutation : nullptr),
        keys_(other.keys_)
  {
    GUDHI_CHECK(!other.is_initialized() || permutation == *other.newToOldPerm_,
                "Only the address of the permutation vector is allowed to change, not its content.");
  }

  /**
   * @brief Copy constructor. Updates the pointers with the given complex and map if @p other was initialized.
   *
   * @param other Interface to copy.
   * @param boundaries Complex. Has to correspond to the same complex than in @p other, except that its address
   * does not have to be the same.
   * @param permutation Permutation map. Has to correspond to the same map than in @p other, except that its address
   * does not have to be the same.
   */
  Multi_parameter_filtered_complex_pcoh_interface(const Multi_parameter_filtered_complex_pcoh_interface &other,
                                                  const Complex &boundaries,
                                                  const Map &permutation)
      : boundaries_(other.is_initialized() ? &boundaries : nullptr),
        newToOldPerm_(other.is_initialized() ? &permutation : nullptr),
        keys_(other.keys_)
  {
    GUDHI_CHECK(!other.is_initialized() || permutation == *other.newToOldPerm_,
                "Only the address of the permutation vector is allowed to change, not its content.");
    GUDHI_CHECK(!other.is_initialized() || boundaries.get_boundaries() == *other.boundaries_->get_boundaries(),
                "Only the address of the complex is allowed to change, not its content.");
  }

  Multi_parameter_filtered_complex_pcoh_interface(Multi_parameter_filtered_complex_pcoh_interface &&other) = delete;

  /**
   * @brief Move constructor. Moves the complex pointer from @p other and updates the permutation pointer with the
   * given map if the complex was initialized.
   *
   * @param other Interface to move.
   * @param permutation Permutation map. Has to correspond to the same map than in @p other, except that its address
   * does not have to be the same.
   */
  Multi_parameter_filtered_complex_pcoh_interface(Multi_parameter_filtered_complex_pcoh_interface &&other,
                                                  const Map &permutation)
      : boundaries_(std::exchange(other.boundaries_, nullptr)),
        newToOldPerm_(boundaries_ != nullptr ? &permutation : nullptr),
        keys_(std::move(other.keys_))
  {
    other.newToOldPerm_ = nullptr;
  }

  /**
   * @brief Move constructor. Updates the pointers with the given complex and map if @p other was initialized.
   *
   * @param other Interface to move.
   * @param boundaries Complex. Has to correspond to the same complex than in @p other, except that its address
   * does not have to be the same.
   * @param permutation Permutation map. Has to correspond to the same map than in @p other, except that its address
   * does not have to be the same.
   */
  Multi_parameter_filtered_complex_pcoh_interface(Multi_parameter_filtered_complex_pcoh_interface &&other,
                                                  const Complex &boundaries,
                                                  const Map &permutation)
      : boundaries_(other.is_initialized() ? &boundaries : nullptr),
        newToOldPerm_(other.is_initialized() ? &permutation : nullptr),
        keys_(std::move(other.keys_))
  {
    other.boundaries_ = nullptr;
    other.newToOldPerm_ = nullptr;
  }

  ~Multi_parameter_filtered_complex_pcoh_interface() = default;

  Multi_parameter_filtered_complex_pcoh_interface &operator=(
      const Multi_parameter_filtered_complex_pcoh_interface &other) = delete;
  Multi_parameter_filtered_complex_pcoh_interface &operator=(
      Multi_parameter_filtered_complex_pcoh_interface &&other) noexcept = delete;

  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_parameter_filtered_complex_pcoh_interface &be1,
                   Multi_parameter_filtered_complex_pcoh_interface &be2) noexcept
  {
    std::swap(be1.boundaries_, be2.boundaries_);
    std::swap(be1.newToOldPerm_, be2.newToOldPerm_);
    be1.keys_.swap(be2.keys_);
  }

  /**
   * @brief Reinitializes the interface with the new given complex and permutation.
   * To use instead of the classical assign operator `operator=`.
   */
  void reinitialize(const Complex &boundaries, const Map &permutation)
  {
    boundaries_ = &boundaries;
    newToOldPerm_ = &permutation;
    keys_ = Map(boundaries.get_number_of_cycle_generators(), -1);
  }

  void reset()
  {
    boundaries_ = nullptr;
    newToOldPerm_ = nullptr;
    keys_.clear();
  }

  /**
   * @brief Returns `true` if and only if all pointers are not null.
   */
  [[nodiscard]] bool is_initialized() const { return boundaries_ != nullptr && newToOldPerm_ != nullptr; }

  [[nodiscard]] std::size_t num_simplices() const { return newToOldPerm_->size(); }

  [[nodiscard]] Filtration_value filtration(Simplex_handle sh) const
  {
    return sh == null_simplex() ? std::numeric_limits<Filtration_value>::max() : keys_[sh];
  }

  [[nodiscard]] Dimension dimension() const { return boundaries_->get_max_dimension(); }

  [[nodiscard]] Dimension dimension(Simplex_handle sh) const
  {
    return sh == null_simplex() ? -1 : boundaries_->get_dimensions()[sh];
  }

  // assumes that pcoh will assign the keys from 0 to n in order of filtrations
  void assign_key(Simplex_handle sh, Simplex_key key)
  {
    if (sh != null_simplex()) keys_[sh] = key;
  }

  [[nodiscard]] Simplex_key key(Simplex_handle sh) const { return sh == null_simplex() ? null_key() : keys_[sh]; }

  static constexpr Simplex_key null_key() { return static_cast<Simplex_key>(-1); }

  [[nodiscard]] Simplex_handle simplex(Simplex_key key) const
  {
    return key == null_key() ? null_simplex() : (*newToOldPerm_)[key];
  }

  static constexpr Simplex_handle null_simplex() { return static_cast<Simplex_handle>(-1); }

  // only used in update_cohomology_groups_edge, so not used without optimizations
  [[nodiscard]] std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) const
  {
    if (sh == null_simplex()) return {null_simplex(), null_simplex()};
    GUDHI_CHECK(dimension(sh) == 1, "Endpoints only available for edges.");
    const auto &col = boundary_simplex_range(sh);
    GUDHI_CHECK(col.size() == 2, "Edge should have two vertices as border.");
    return {col[0], col[1]};
  }

  [[nodiscard]] const Filtration_simplex_range &filtration_simplex_range() const { return *newToOldPerm_; }

  [[nodiscard]] const Boundary_simplex_range &boundary_simplex_range(Simplex_handle sh) const
  {
    return boundaries_->get_boundaries()[sh];
  }

  friend std::ostream &operator<<(std::ostream &stream,
                                  const Multi_parameter_filtered_complex_pcoh_interface &complex)
  {
    stream << "[\n";
    for (auto i : complex.filtration_simplex_range()) {
      stream << "[";
      for (const auto &idx : complex.boundary_simplex_range(i)) stream << complex.keys_[idx] << ", ";
      stream << "]\n";
    }

    stream << "]\n";
    return stream;
  }

 private:
  Complex const *boundaries_; /**< Pointer to complex. */
  Map const *newToOldPerm_;   /**< Pointer to filtration position to complex position map. */
  Map keys_;                  /**< Keys assigned to a cell. TODO: potentially the identity. If yes, to remove. */
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_COMPLEX_PCOH_INTERFACE_H_INCLUDED
