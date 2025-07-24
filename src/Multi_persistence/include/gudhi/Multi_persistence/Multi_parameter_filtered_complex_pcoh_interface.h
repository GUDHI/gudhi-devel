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

// #include <cstddef>
#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_parameter_filtered_complex.h>

namespace Gudhi {
namespace multi_persistence {

template <class MultiFiltrationValue>
class Multi_parameter_filtered_complex_pcoh_interface
{
 public:
  using Complex = Multi_parameter_filtered_complex<MultiFiltrationValue>;
  using Simplex_key = std::uint32_t;
  using Simplex_handle = Simplex_key;
  using Filtration_value = Simplex_key;
  using Dimension = int;
  using Map = std::vector<Simplex_handle>;
  using Filtration_simplex_range = Map;
  using Boundary_simplex_range = Map;

  Multi_parameter_filtered_complex_pcoh_interface() : boundaries_(nullptr), newToOldPerm_(nullptr) {}

  // `boundaries` and `permutation` are assumed to have stable size, i.e., their address never changes
  Multi_parameter_filtered_complex_pcoh_interface(const Complex &boundaries, const Map &permutation)
      : boundaries_(&boundaries),
        newToOldPerm_(&permutation),
        keys_(boundaries.get_number_of_cycle_generators(), -1) /* ,
         keysToSH_(permutation.size(), -1) */
  {}

  Multi_parameter_filtered_complex_pcoh_interface(const Multi_parameter_filtered_complex_pcoh_interface &other) =
      delete;

  // permutation and boundaries is assumed to be the same than from the copied object, just its address can change
  Multi_parameter_filtered_complex_pcoh_interface(const Multi_parameter_filtered_complex_pcoh_interface &other,
                                                  const Map &permutation)
      : boundaries_(other.boundaries_), newToOldPerm_(other.is_initialized() ? &permutation : nullptr), keys_(other.keys_) /* ,
                                                                keysToSH_(other.keysToSH_) */
  {
    GUDHI_CHECK(!other.is_initialized() || permutation == *other.newToOldPerm_,
                "Only the address of the permutation vector is allowed to change, not its content.");
  }

  // permutation and boundaries is assumed to be the same than from the copied object, just its address can change
  Multi_parameter_filtered_complex_pcoh_interface(const Multi_parameter_filtered_complex_pcoh_interface &other,
                                                  const Complex &boundaries,
                                                  const Map &permutation)
      : boundaries_(other.is_initialized() ? &boundaries : nullptr), newToOldPerm_(other.is_initialized() ? &permutation : nullptr), keys_(other.keys_) /* ,
                                                                keysToSH_(other.keysToSH_) */
  {
    GUDHI_CHECK(!other.is_initialized() || permutation == *other.newToOldPerm_,
                "Only the address of the permutation vector is allowed to change, not its content.");
    GUDHI_CHECK(!other.is_initialized() || boundaries.get_boundaries() == *other.boundaries_->get_boundaries(),
                "Only the address of the complex is allowed to change, not its content.");
  }

  Multi_parameter_filtered_complex_pcoh_interface(Multi_parameter_filtered_complex_pcoh_interface &&other) = delete;

  // permutation is assumed to be the same than from the moved object, just its address can change
  Multi_parameter_filtered_complex_pcoh_interface(Multi_parameter_filtered_complex_pcoh_interface &&other,
                                                  const Map &permutation)
      : boundaries_(std::exchange(other.boundaries_, nullptr)),
        newToOldPerm_(other.is_initialized() ? &permutation : nullptr),
        keys_(std::move(other.keys_)) /* ,
keysToSH_(std::move(other.keysToSH_)) */
  {
    other.newToOldPerm_ = nullptr;
  }

  // permutation is assumed to be the same than from the moved object, just its address can change
  Multi_parameter_filtered_complex_pcoh_interface(Multi_parameter_filtered_complex_pcoh_interface &&other,
                                                  const Complex &boundaries,
                                                  const Map &permutation)
      : boundaries_(other.is_initialized() ? &boundaries : nullptr), newToOldPerm_(other.is_initialized() ? &permutation : nullptr), keys_(std::move(other.keys_)) /* ,
                                                                keysToSH_(std::move(other.keysToSH_)) */
  {
    other.boundaries_ = nullptr;
    other.newToOldPerm_ = nullptr;
  }

  friend void swap(Multi_parameter_filtered_complex_pcoh_interface &be1,
                   Multi_parameter_filtered_complex_pcoh_interface &be2)
  {
    std::swap(be1.boundaries_, be2.boundaries_);
    std::swap(be1.newToOldPerm_, be2.newToOldPerm_);
    be1.keys_.swap(be2.keys_);
    // be1.keysToSH_.swap(be2.keysToSH_);
  }

  void reinitialize(const Complex &boundaries, const Map &permutation)
  {
    boundaries_ = &boundaries;
    newToOldPerm_ = &permutation;
    keys_ = Map(boundaries.get_number_of_cycle_generators(), -1);
    // keysToSH_ = Map(permutation.size(), -1);
  }

  bool is_initialized() const { return boundaries_ != nullptr && newToOldPerm_ != nullptr; }

  std::size_t num_simplices() const { return newToOldPerm_->size(); }

  Filtration_value filtration(Simplex_handle sh) const
  {
    return sh == null_simplex() ? std::numeric_limits<Filtration_value>::max() : keys_[sh];
  }

  Dimension dimension() const { return boundaries_->get_max_dimension(); }

  Dimension dimension(Simplex_handle sh) const { return sh == null_simplex() ? -1 : boundaries_->get_dimensions()[sh]; }

  // assumes that pcoh will assign the keys from 0 to n in order of filtrations
  void assign_key(Simplex_handle sh, Simplex_key key)
  {
    if (sh != null_simplex()) keys_[sh] = key;
    // if (key != null_key()) keysToSH_[key] = sh;
  }

  Simplex_key key(Simplex_handle sh) const { return sh == null_simplex() ? null_key() : keys_[sh]; }

  static constexpr Simplex_key null_key() { return static_cast<Simplex_key>(-1); }

  // Simplex_handle simplex(Simplex_key key) const { return key == null_key() ? null_simplex() : keysToSH_[key]; }
  Simplex_handle simplex(Simplex_key key) const { return key == null_key() ? null_simplex() : (*newToOldPerm_)[key]; }

  static constexpr Simplex_handle null_simplex() { return static_cast<Simplex_handle>(-1); }

  // only used in update_cohomology_groups_edge, so not used without optimizations
  std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) const
  {
    if (sh == null_simplex()) return {null_simplex(), null_simplex()};
    GUDHI_CHECK(dimension(sh) == 1, "Endpoints only avalaible for edges.");
    const auto &col = boundary_simplex_range(sh);
    GUDHI_CHECK(col.size() == 2, "Edge should have two vertices as border.");
    return {col[0], col[1]};
  }

  const Filtration_simplex_range &filtration_simplex_range() const { return *newToOldPerm_; }

  const Boundary_simplex_range &boundary_simplex_range(Simplex_handle sh) const
  {
    return boundaries_->get_boundaries()[sh];
  }

  friend std::ostream &operator<<(std::ostream &stream,
                                  const Multi_parameter_filtered_complex_pcoh_interface &interface)
  {
    // std::vector<Simplex_key> inv(interface.newToOldPerm_->size());
    // for (unsigned int i = 0; i < interface.newToOldPerm_->size(); ++i) {
    //   inv[(*interface.newToOldPerm_)[i]] = i;
    // }

    stream << "[\n";
    for (auto i : interface.filtration_simplex_range()) {
      stream << "[";
      // for (const auto &stuff : interface.boundary_simplex_range(i)) stream << inv[stuff] << ", ";
      for (const auto &stuff : interface.boundary_simplex_range(i)) stream << interface.keys_[stuff] << ", ";
      stream << "]\n";
    }

    stream << "]\n";
    return stream;
  }

 private:
  Complex const *boundaries_;
  Map const *newToOldPerm_;
  Map keys_;
  // Map keysToSH_;
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_COMPLEX_PCOH_INTERFACE_H_INCLUDED
