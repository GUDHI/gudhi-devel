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

#include <ostream>  //std::ostream
#include <limits>   //std::numeric_limits
#include <tuple>
#include <utility>  //std::move

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Persistence_interval persistence_interval.h gudhi/persistence_interval.h
 * @ingroup persistence_matrix
 *
 * @brief Type for an interval in a persistent diagram or barcode.
 * Stores the birth, death and dimension of the interval.
 * 
 * @tparam dimension_type Type of the dimension value.
 * @tparam event_value_type Type of the birth and death value.
 */
template <typename dimension_type, typename event_value_type>
struct Persistence_interval {
  /**
   * @brief Stores the infinity value for birth and death events. Its value depends on the template parameter
   * `event_value_type`:
   * - if `event_value_type` has a native infinity value, it takes this value,
   * - otherwise, if `event_value_type` is a signed type, it takes value -1,
   * - otherwise, if `event_value_type` is a unsigned type, it takes the maximal possible value.
   *
   * Is also used as default value for birth and death attributes when not initialized.
   */
  static constexpr event_value_type inf = std::numeric_limits<event_value_type>::has_infinity
                                              ? std::numeric_limits<event_value_type>::infinity()
                                              : static_cast<event_value_type>(-1);

  /**
   * @brief Default constructor. Initializes the stored dimension to -1 and the stored birth and death values
   * to @ref inf.
   */
  Persistence_interval() : dim(-1), birth(inf), death(inf) {}

  /**
   * @brief Constructor. Initializes the stored dimension and the stored birth to the given values and the stored
   * death value to @ref inf.
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

/**
 * @ingroup persistence_matrix
 *
 * @brief Overload of `std::tuple_size` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam dimension_type First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam event_value_type Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 */
template <typename dimension_type, typename event_value_type>
struct tuple_size<Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type> >
    : std::integral_constant<std::size_t, 3> {};

/**
 * @ingroup persistence_matrix
 *
 * @brief Overload of `std::tuple_element` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam I Index of the type to store: 0 for the birth value type, 1 for the death value type and 2 for the
 * dimension value type.
 * @tparam dimension_type First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam event_value_type Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 */
template <std::size_t I, typename dimension_type, typename event_value_type>
struct tuple_element<I, Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type> > {
  static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

  using type = typename std::conditional<I < 2, event_value_type, dimension_type>::type;
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Overload of `std::get` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
 * @tparam dimension_type First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam event_value_type Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @param i Interval from which the value should be returned.
 * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
 */
template <size_t I, typename dimension_type, typename event_value_type>
constexpr auto& get(Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>& i) noexcept {
  static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

  if constexpr (I == 0) return i.birth;
  if constexpr (I == 1) return i.death;
  if constexpr (I == 2) return i.dim;
}

/**
 * @ingroup persistence_matrix
 *
 * @brief Overload of `std::get` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
 * @tparam dimension_type First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam event_value_type Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @param i Interval from which the value should be returned.
 * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
 */
template <size_t I, typename dimension_type, typename event_value_type>
constexpr const auto& get(
    const Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>& i) noexcept {
  static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

  if constexpr (I == 0) return i.birth;
  if constexpr (I == 1) return i.death;
  if constexpr (I == 2) return i.dim;
}

/**
 * @ingroup persistence_matrix
 *
 * @brief Overload of `std::get` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
 * @tparam dimension_type First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam event_value_type Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @param i Interval from which the value should be returned.
 * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
 */
template <size_t I, typename dimension_type, typename event_value_type>
constexpr auto&& get(Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>&& i) noexcept {
  static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

  if constexpr (I == 0) return std::move(i.birth);
  if constexpr (I == 1) return std::move(i.death);
  if constexpr (I == 2) return std::move(i.dim);
}

/**
 * @ingroup persistence_matrix
 *
 * @brief Overload of `std::get` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
 * @tparam dimension_type First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam event_value_type Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @param i Interval from which the value should be returned.
 * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
 */
template <size_t I, typename dimension_type, typename event_value_type>
constexpr const auto&& get(
    const Gudhi::persistence_matrix::Persistence_interval<dimension_type, event_value_type>&& i) noexcept {
  static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

  if constexpr (I == 0) return std::move(i.birth);
  if constexpr (I == 1) return std::move(i.death);
  if constexpr (I == 2) return std::move(i.dim);
}

}  // namespace std

#endif  // PM_INTERVAL_INCLUDED
