/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/04 Hannah Schreiber: Reorganization.
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Thread_safe_slicer.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Thread_safe_slicer class.
 */

#ifndef MP_THREAD_SAFE_SLICER_H_INCLUDED
#define MP_THREAD_SAFE_SLICER_H_INCLUDED

#include <ostream>
#include <vector>

namespace Gudhi {
namespace multi_persistence {

template <class Slicer>
class Thread_safe_slicer : private Slicer
{
 public:
  using Persistence = typename Slicer::Persistence;
  using Filtration_value = typename Slicer::Filtration_value;
  using T = typename Slicer::T;
  using Complex = typename Slicer::Complex;
  using Index = typename Slicer::Index;
  using Dimension = typename Slicer::Dimension;
  template <typename Value = T>
  using Bar = typename Slicer::template Bar<Value>;
  template <typename Value = T>
  using Barcode = typename Slicer::template Barcode<Value>;
  template <typename Value = T>
  using Flat_barcode = typename Slicer::template Flat_barcode<Value>;
  template <typename Value = T>
  using Multi_dimensional_barcode = typename Slicer::template Multi_dimensional_barcode<Value>;
  template <typename Value = T>
  using Multi_dimensional_flat_barcode = typename Slicer::template Multi_dimensional_flat_barcode<Value>;
  using Cycle = typename Slicer::Cycle;
  using Thread_safe = Thread_safe_slicer;

  // CONSTRUCTORS

  Thread_safe_slicer(const Slicer& slicer)
      : Slicer(slicer.get_slice(), slicer.get_current_order(), slicer.get_persistence_algorithm()), slicer_(&slicer)
  {}

  Thread_safe_slicer(const Thread_safe_slicer& slicer)
      : Slicer(slicer.get_slice(), slicer.get_current_order(), slicer.get_persistence_algorithm()),
        slicer_(slicer.slicer_)
  {}

  // ACCESS

  Thread_safe weak_copy() const { return Thread_safe_slicer(*this); }

  Index get_number_of_cycle_generators() const { return slicer_->get_number_of_cycle_generators(); }

  Index get_number_of_parameters() const { return slicer_->get_number_of_parameters(); }

  const std::vector<Index>& get_current_order() const { return Slicer::get_current_order(); }

  const std::vector<T>& get_slice() const { return Slicer::get_slice(); }

  const Persistence& get_persistence_algorithm() const { return Slicer::get_persistence_algorithm(); }

  std::pair<Filtration_value, Filtration_value> get_bounding_box() const { return slicer_->get_bounding_box(); }

  const typename Complex::Filtration_value_container& get_filtration_values() const
  {
    return slicer_->get_filtration_values();
  }

  const std::vector<Dimension>& get_dimensions() const { return slicer_->get_dimensions(); }

  const typename Complex::Boundary_container& get_boundaries() const { return slicer_->get_boundaries(); }

  Dimension get_dimension(Index i) const { return slicer_->get_dimension(i); }

  // MODIFIERS

  template <class Array = std::initializer_list<T>>
  void set_slice(const Array& slice)
  {
    Slicer::set_slice(slice);
  }

  template <class Line>
  void push_to(const Line& line)
  {
    Slicer::_push_to(slicer_->complex_, line);
  }

  // PERSISTENCE

  bool persistence_computation_is_initialized() const { return Slicer::persistence_computation_is_initialized(); }

  void initialize_persistence_computation(const bool ignoreInf = true)
  {
    Slicer::_initialize_persistence_computation(slicer_->complex_, ignoreInf);
  }

  void vineyard_update() { Slicer::vineyard_update(); }

  template <typename Value = T>
  auto get_barcode(int maxDim = -1)
  {
    // complex in parent is empty, so maxDim needs to be initialized from the outside.
    if (maxDim < 0) maxDim = slicer_->get_max_dimension();
    return Slicer::get_barcode(maxDim);
  }

  template <bool withDim = false, typename Value = T>
  auto get_flat_barcode(int maxDim = -1)
  {
    // complex in parent is empty, so maxDim needs to be initialized from the outside.
    if (maxDim < 0) maxDim = slicer_->get_max_dimension();
    return Slicer::get_flat_barcode(maxDim);
  }

  std::vector<Multi_dimensional_barcode<T>> persistence_on_lines(const std::vector<std::vector<T>>& basePoints,
                                                                 bool ignoreInf)
  {
    return Slicer::persistence_on_lines(basePoints, ignoreInf);
  }

  std::vector<Multi_dimensional_barcode<T>> persistence_on_lines(
      const std::vector<std::pair<std::vector<T>, std::vector<T>>>& basePointsWithDirections,
      bool ignoreInf)
  {
    return Slicer::persistence_on_lines(basePointsWithDirections, ignoreInf);
  }

  std::vector<std::vector<Cycle>> get_representative_cycles(bool update = true)
  {
    return Slicer::_get_representative_cycles(slicer_->complex_, update);
  }

  // FRIENDS

  friend std::ostream& operator<<(std::ostream& stream, const Thread_safe_slicer& slicer)
  {
    stream << "-------------------- Thread_safe_slicer \n";

    stream << "--- Filtered complex \n";
    stream << slicer.slicer_->complex_;

    stream << "--- Order \n";
    stream << "{";
    for (const auto& idx : slicer.get_current_order()) stream << idx << ", ";
    stream << "}" << std::endl;

    stream << "--- Current slice filtration\n";
    stream << "{";
    for (const auto& val : slicer.get_slice()) stream << val << ", ";
    stream << "\b" << "\b";
    stream << "}" << std::endl;

    stream << "--- PersBackend \n";
    stream << slicer.get_persistence_algorithm();

    return stream;
  }

 private:
  Slicer const* slicer_;
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_THREAD_SAFE_SLICER_H_INCLUDED
