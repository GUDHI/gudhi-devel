/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/04 Hannah Schreiber: Reorganization + documentation.
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
#include <utility>
#include <vector>

#include <gudhi/Multi_persistence/Line.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Thread_safe_slicer Thread_safe_slicer.h gudhi/Thread_safe_slicer.h
 * @ingroup multi_persistence
 *
 * @brief A more or less "thread safe version" of @ref Slicer. Gives access to all its const methods and persistence
 * related methods. Each instance will store a pointer to the original slicer, but will have its own copies of all
 * persistence related containers. It corresponds therefore more to a "light/week copy" of a slicer.
 *
 * @tparam Slicer Underlying @ref Slicer type.
 */
template <class Slicer>
class Thread_safe_slicer : private Slicer
{
 public:
  using Persistence = typename Slicer::Persistence;           /**< Persistence algorithm type. */
  using Filtration_value = typename Slicer::Filtration_value; /**< Filtration value type. */
  using T = typename Slicer::T;                               /**< Numerical filtration value element type. */
  using Complex = typename Slicer::Complex;                   /**< Complex type. */
  using Index = typename Slicer::Index;                       /**< Complex index type. */
  using Dimension = typename Slicer::Dimension;               /**< Dimension type. */
  template <typename Value = T>
  using Bar = typename Slicer::template Bar<Value>;           /**< Bar type. */
  /**
   * @brief Barcode type. A vector of @ref Bar, a tuple like structure containing birth, death and dimension of a bar.
   */
  template <typename Value = T>
  using Barcode = typename Slicer::template Barcode<Value>;
  /**
   * @brief Flat barcode type. All bars are represented by a birth and a death value stored respectively at even and
   * odd indices of the vector.
   */
  template <typename Value = T>
  using Flat_barcode = typename Slicer::template Flat_barcode<Value>;
  /**
   * @brief Barcode ordered by dimension type. A vector which has at index \f$ d \f$ the @ref Barcode of dimension
   * \f$ d \f$.
   */
  template <typename Value = T>
  using Multi_dimensional_barcode = typename Slicer::template Multi_dimensional_barcode<Value>;
  /**
   * @brief Flat barcode ordered by dimension type. A vector which has at index \f$ d \f$ the @ref Flat_barcode of
   * dimension \f$ d \f$.
   */
  template <typename Value = T>
  using Multi_dimensional_flat_barcode = typename Slicer::template Multi_dimensional_flat_barcode<Value>;
  using Cycle = typename Slicer::Cycle;   /**< Cycle type. */
  using Thread_safe = Thread_safe_slicer; /**< Thread safe slicer type. */

  // CONSTRUCTORS

  /**
   * @brief Constructor. Will store a pointer to the given slicer and copy all persistence related container.
   *
   * @param slicer Original slicer.
   */
  Thread_safe_slicer(const Slicer& slicer)
      : Slicer(slicer.get_slice(), slicer.get_current_order(), slicer.get_persistence_algorithm()), slicer_(&slicer)
  {}

  /**
   * @brief Copy constructor.
   */
  Thread_safe_slicer(const Thread_safe_slicer& slicer)
      : Slicer(slicer.get_slice(), slicer.get_current_order(), slicer.get_persistence_algorithm()),
        slicer_(slicer.slicer_)
  {}

  /**
   * @brief Move constructor.
   */
  Thread_safe_slicer(Thread_safe_slicer&& slicer) noexcept
      : Slicer(std::move(slicer.slice_), std::move(slicer.generatorOrder_), std::move(slicer.persistence_)),
        slicer_(slicer.slicer_)
  {}

  ~Thread_safe_slicer() = default;

  /**
   * @brief Assign operator.
   */
  Thread_safe_slicer& operator=(const Thread_safe_slicer& other)
  {
    if (this == &other) return *this;

    Slicer::slice_ = other.slice_;
    Slicer::generatorOrder_ = other.generatorOrder_;
    Slicer::persistence_.reinitialize(slicer_->complex_, Slicer::generatorOrder_);
    slicer_ = other.slicer_;

    return *this;
  }

  /**
   * @brief Move assign operator.
   */
  Thread_safe_slicer& operator=(Thread_safe_slicer&& other) noexcept = delete;

  // ACCESS

  /**
   * @brief Returns a copy of this object.
   */
  Thread_safe weak_copy() const { return Thread_safe_slicer(*this); }

  /**
   * @brief Returns the number of generators in the stored module.
   */
  Index get_number_of_cycle_generators() const { return slicer_->get_number_of_cycle_generators(); }

  /**
   * @brief Returns the number of parameters of the stored filtration. If the module is empty, the number returned is 0.
   */
  Index get_number_of_parameters() const { return slicer_->get_number_of_parameters(); }

  /**
   * @brief Returns a const reference to the current permutation map, indicating in which order are the generators
   * with respect to the current slice (i.e., \$f order[i] \$f corresponds to the index in the complex of the
   * \$f i^{th} \$f generator in the filtration represented by the slice). It will be initialized with
   * @ref initialize_persistence_computation.
   *
   * If `ignoreInf` was true when calling @ref initialize_persistence_computation, indices of generators at infinity
   * are not stored in the container. That means that the size can be smaller than what
   * @ref get_number_of_cycle_generators returns.
   */
  const std::vector<Index>& get_current_order() const { return Slicer::get_current_order(); }

  /**
   * @brief Returns a const reference to the current slice. It can be initialized or updated with @ref set_slice
   * and @ref push_to.
   */
  const std::vector<T>& get_slice() const { return Slicer::get_slice(); }

  /**
   * @brief Returns a reference to the current slice. It can also be initialized or updated with @ref set_slice
   * and @ref push_to.
   */
  std::vector<T>& get_slice() { return Slicer::get_slice(); }

  /**
   * @brief Returns a const reference to the class computing the persistence of the current slice. It will be
   * initialized with @ref initialize_persistence_computation.
   */
  const Persistence& get_persistence_algorithm() const { return Slicer::get_persistence_algorithm(); }

  /**
   * @brief Returns two filtration values representing respectively the greatest common lower bound of all filtration
   * values in the filtration and the lowest common upper bound of them.
   */
  std::pair<Filtration_value, Filtration_value> get_bounding_box() const { return slicer_->get_bounding_box(); }

  /**
   * @brief Returns a const reference to the filtration value container. A filtration value at index \$f i \$f
   * correspond to the filtration value associated to the generators at index \$f i \$f.
   */
  const typename Complex::Filtration_value_container& get_filtration_values() const
  {
    return slicer_->get_filtration_values();
  }

  /**
   * @brief Returns a const reference to the filtration value associated to the generator at index \$f i \$f.
   */
  const Filtration_value& get_filtration_value(Index i) const { return slicer_->get_filtration_value(i); }

  /**
   * @brief Returns a const reference to the dimension container. A value at index \$f i \$f corresponds to the
   * dimension of the generator at index \$f i \$f.
   */
  const std::vector<Dimension>& get_dimensions() const { return slicer_->get_dimensions(); }

  /**
   * @brief Returns the dimension of the generator at index \$f i \$f.
   */
  Dimension get_dimension(Index i) const { return slicer_->get_dimension(i); }

  /**
   * @brief Returns the maximal dimension of a generator in the module.
   */
  Dimension get_max_dimension() const { return slicer_->get_max_dimension(); }

  /**
   * @brief Returns a const reference to the boundary container. The element at index \$f i \$f corresponds to the
   * boundary of the generator at index \$f i \$f.
   */
  const typename Complex::Boundary_container& get_boundaries() const { return slicer_->get_boundaries(); }

  /**
   * @brief Returns the boundary of the generator at index \$f i \$f.
   */
  const typename Complex::Boundary& get_boundary(Index i) const { return slicer_->get_boundary(i); }

  // MODIFIERS

  /**
   * @brief Sets the current slice, that is the 1-parameter filtration values associated to each generator on that line.
   * The value at \$f slice[i] \$f has to corresponds to the value for the generator at index \$f i \$f.
   * One can also sets the slice directly from the line with @ref push_to.
   *
   * @tparam Array Container which can be converted into a vector of `T`.
   */
  template <class Array = std::initializer_list<T>>
  void set_slice(const Array& slice)
  {
    Slicer::set_slice(slice);
  }

  /**
   * @brief Sets the current slice by computing the 1-parameter filtration values fo each generator on the given line.
   *
   * @tparam Line_like Any type convertible to a @ref Line class. Default value: `std::initializer_list<T>`.
   */
  template <class Line_like = std::initializer_list<T>>
  void push_to(const Line_like& line)
  {
    Slicer::_push_to(slicer_->complex_, Line<typename Line_like::value_type>(line));
  }

  /**
   * @brief Sets the current slice by computing the 1-parameter filtration values fo each generator on the given line.
   *
   * @tparam U Template parameter of the given line.
   */
  template <class U>
  void push_to(const Line<U>& line)
  {
    Slicer::_push_to(slicer_->complex_, line);
  }

  // PERSISTENCE

  /**
   * @brief Returns true if and only if @ref initialize_persistence_computation was properly called.
   */
  [[nodiscard]] bool persistence_computation_is_initialized() const
  {
    return Slicer::persistence_computation_is_initialized();
  }

  /**
   * @brief Initializes the persistence computation of the current slice. If the slice was not set properly as
   * a valid 1-dimensional filtration, the behaviour is undefined.
   *
   * @param ignoreInf If true, all cells at infinity filtration values are ignored for the initialization, resulting
   * potentially in less storage use and better performance. But note that this can be problematic with the use of
   * @ref vineyard_update. Default value: true.
   */
  void initialize_persistence_computation(const bool ignoreInf = true)
  {
    Slicer::_initialize_persistence_computation(slicer_->complex_, ignoreInf);
  }

  /**
   * @brief After the persistence computation was initialized for a slice and the slice changes, this method can
   * update everything necessary for the barcode without re-computing everything from scratch (contrary to
   * @ref initialize_persistence_computation). Furthermore, it guarantees that the new barcode will "match" the
   * precedent one. TODO: explain exactly what it means and how to do the matching.
   * The method will have better performance if the complex is ordered by dimension.
   *
   * Only available if PersistenceAlgorithm::is_vine is true.
   *
   * @pre @ref initialize_persistence_computation has to be called at least once before.
   *
   * @warning If `ignoreInf` was set to true when initializing the persistence computation, any update of the slice has
   * to keep at infinity the boundaries which were before, otherwise the behaviour is undefined (it will throw with
   * high probability).
   */
  void vineyard_update() { Slicer::vineyard_update(); }

  /**
   * @brief Returns the barcode of the current slice. The barcode format will change depending on the template values.
   *
   * @pre @ref initialize_persistence_computation has to be called at some point before.
   *
   * @tparam byDim If true, the barcode is returned as @ref Multi_dimensional_barcode, otherwise as @ref Barcode.
   * @tparam Value Type of the birth and death values.
   * @param maxDim Maximal dimension to be included in the barcode. If negative, all dimensions are included.
   * Default value: -1.
   */
  template <bool byDim = true, typename Value = T, bool idx = false>
  auto get_barcode(int maxDim = -1)
  {
    // complex in parent is empty, so maxDim needs to be initialized from the outside.
    if (maxDim < 0) maxDim = slicer_->get_max_dimension();
    return Slicer::template get_barcode<byDim, Value, idx>(maxDim);
  }

  /**
   * @brief Returns the barcode of the current slice. The barcode format will change depending on the template values.
   *
   * @pre @ref initialize_persistence_computation has to be called at some point before.
   *
   * @tparam byDim If true, the barcode is returned as @ref Multi_dimensional_flat_barcode, otherwise as
   * @ref Flat_barcode.
   * @tparam Value Type of the birth and death values.
   * @param maxDim Maximal dimension to be included in the barcode. If negative, all dimensions are included.
   * Default value: -1.
   */
  template <bool byDim = false, typename Value = T, bool idx = false>
  auto get_flat_barcode(int maxDim = -1)
  {
    // complex in parent is empty, so maxDim needs to be initialized from the outside.
    if (maxDim < 0) maxDim = slicer_->get_max_dimension();
    return Slicer::template get_flat_barcode<byDim, Value, idx>(maxDim);
  }

  /**
   * @brief Returns the representative cycles of the current slice. All cycles of dimension \f$ d \f$ are stored at
   * index \f$ d \f$ of the returned vector. A cycle is represented by a vector of boundary indices. That is, the index
   * \f$ i \f$ in a cycle represents the cell which boundary can be retrieved by @ref get_boundary "get_boundary(i)".
   *
   * Only available if PersistenceAlgorithm::has_rep_cycles is true.
   *
   * @pre @ref initialize_persistence_computation has to be called at least once before.
   *
   * @param update If true, updates the stored representative cycles, otherwise just returns the container in its
   * current state. So should be true at least the first time the method is used.
   */
  std::vector<std::vector<Cycle>> get_representative_cycles(bool update = true)
  {
    return Slicer::_get_representative_cycles(slicer_->complex_, update);
  }

  Cycle get_most_persistent_cycle(Dimension dim = 1, bool update = true)
  {
    return Slicer::get_most_persistent_cycle(dim, update);
  }

  // FRIENDS

  /**
   * @brief Outstream operator.
   */
  friend std::ostream& operator<<(std::ostream& stream, Thread_safe_slicer& slicer)
  {
    stream << "-------------------- Thread_safe_slicer \n";

    stream << "--- Filtered complex \n";
    stream << slicer.slicer_->complex_;

    stream << "--- Order \n";
    stream << "{";
    for (const auto& idx : slicer.get_current_order()) stream << idx << ", ";
    stream << "}" << '\n';

    stream << "--- Current slice filtration\n";
    stream << "{";
    for (const auto& val : slicer.get_slice()) stream << val << ", ";
    stream << "\b" << "\b";
    stream << "}" << '\n';

    stream << "--- PersBackend \n";
    stream << slicer.persistence_;

    return stream;
  }

 private:
  Slicer const* slicer_; /**< Original slicer. */
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_THREAD_SAFE_SLICER_H_INCLUDED
