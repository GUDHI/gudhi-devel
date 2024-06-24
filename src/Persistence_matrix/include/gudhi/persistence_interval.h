/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file persistence_interval.h
 * @author Hannah Schreiber
 * @brief Contains @ref Gudhi::persistence_matrix::Persistence_interval class.
 */

#ifndef PM_INTERVAL_INCLUDED
#define PM_INTERVAL_INCLUDED

#include <ostream>
#include <limits>
#include <tuple>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Type for an interval in a persistent diagram or barcode.
 * Stores the birth, death and dimension of the interval.
 */
template <typename dimension_type, typename event_value_type>
struct Persistence_interval {
  /**
   * @brief Stores the infinity value for birth and death events.
   *
   * Is also used as default value for birth and death attributes when not initialized.
   */
  static constexpr event_value_type inf = std::numeric_limits<event_value_type>::has_infinity
                                              ? std::numeric_limits<event_value_type>::infinity()
                                              : static_cast<event_value_type>(-1);

  /**
   * @brief Default constructor. Initializes the stored dimension to -1 and the stored birth and death values
   * to @ref Persistence_interval::inf.
   */
  Persistence_interval() : dim(-1), birth(inf), death(inf) {}

  /**
   * @brief Constructor. Initializes the stored dimension and the stored birth to the given values and the stored
   * death value to @ref Persistence_interval::inf.
   * 
   * @param dim Dimension of the cycle.
   * @param birth Birth value of the cycle.
   */
  Persistence_interval(dimension_type dim, event_value_type birth) : dim(dim), birth(birth), death(inf) {}

  /**
   * @brief Constructor. Initializes the stored dimension, the stored birth value and the stored
   * death value to the given values.
   * 
   * @param dim Dimension of the cycle.
   * @param birth Birth value of the cycle.
   * @param death Death value of the cycle.
   */
  Persistence_interval(dimension_type dim, event_value_type birth, event_value_type death)
      : dim(dim), birth(birth), death(death) {}

  dimension_type dim;     /**< Dimension of the cycle.*/
  event_value_type birth; /**< Birth value of the cycle. */
  event_value_type death; /**< Death value of the cycle. */

  inline friend std::ostream &operator<<(std::ostream &stream, const Persistence_interval &interval) {
    stream << "[" << interval.dim << "] ";
    stream << interval.birth << ", " << interval.death;
    return stream;
  }
};

}  // namespace persistence_matrix
}  // namespace Gudhi

namespace std {

template <typename dimension_type, typename event_value_type>
struct tuple_size<Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type> >
    : std::integral_constant<std::size_t, 3> {};

template <std::size_t I, typename dimension_type, typename event_value_type>
struct tuple_element<I, Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type> > {
  static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

  using type = typename std::conditional<I < 2, event_value_type, dimension_type>::type;
};

template <size_t I, typename dimension_type, typename event_value_type>
constexpr typename tuple_element<I, tuple<event_value_type, event_value_type, dimension_type> >::type& get(
    Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>& i) noexcept {
  if constexpr (I == 0) return i.birth;
  if constexpr (I == 1) return i.death;
  if constexpr (I == 2) return i.dim;
  // does not compile if I > 2, because of tuple_element
}

template <size_t I, typename dimension_type, typename event_value_type>
constexpr const typename tuple_element<I, tuple<event_value_type, event_value_type, dimension_type> >::type& get(
    const Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>& i) noexcept {
  if constexpr (I == 0) return i.birth;
  if constexpr (I == 1) return i.death;
  if constexpr (I == 2) return i.dim;
  // does not compile if I > 2, because of tuple_element
}

template <size_t I, typename dimension_type, typename event_value_type>
constexpr typename tuple_element<I, tuple<event_value_type, event_value_type, dimension_type> >::type&& get(
    Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>&& i) noexcept {
  if constexpr (I == 0) return i.birth;
  if constexpr (I == 1) return i.death;
  if constexpr (I == 2) return i.dim;
  // does not compile if I > 2, because of tuple_element
}

template <size_t I, typename dimension_type, typename event_value_type>
constexpr const typename tuple_element<I, tuple<event_value_type, event_value_type, dimension_type> >::type&& get(
    const Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>&& i) noexcept {
  if constexpr (I == 0) return i.birth;
  if constexpr (I == 1) return i.death;
  if constexpr (I == 2) return i.dim;
  // does not compile if I > 2, because of tuple_element
}

template <typename dimension_type, typename event_value_type>
constexpr dimension_type& get(
    Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>& t) noexcept {
  return t.dim;
}

template <typename dimension_type, typename event_value_type>
constexpr dimension_type&& get(
    Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>&& t) noexcept {
  return t.dim;
}

template <typename dimension_type, typename event_value_type>
constexpr const dimension_type& get(
    const Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>& t) noexcept {
  return t.dim;
}

template <typename dimension_type, typename event_value_type>
constexpr const dimension_type&& get(
    const Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>&& t) noexcept {
  return t.dim;
}

}  // namespace std

#endif  // PM_INTERVAL_INCLUDED
