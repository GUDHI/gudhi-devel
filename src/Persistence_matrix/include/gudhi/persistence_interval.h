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
 * @brief Type for an interval in a persistent diagram or barcode. Stores the birth, death and dimension of the
 * interval. It can be used as a tuple with `get`/`std::get` (birth, death and dimension in this order),
 * `std::tuple_element` and `std::tuple_size`, as well as structured binding.
 *
 * @tparam Dimension Type of the dimension value.
 * @tparam Event_value Type of the birth and death value.
 */
template <typename Dimension, typename Event_value>
struct Persistence_interval {
  /**
   * @brief Stores the infinity value for birth and death events. Its value depends on the template parameter
   * `Event_value`:
   * - if `Event_value` has a native infinity value, it takes this value,
   * - otherwise, if `Event_value` is a signed type, it takes value -1,
   * - otherwise, if `Event_value` is a unsigned type, it takes the maximal possible value.
   *
   * Is also used as default value for birth and death attributes when not initialized.
   */
  static constexpr Event_value inf = std::numeric_limits<Event_value>::has_infinity
                                         ? std::numeric_limits<Event_value>::infinity()
                                         : static_cast<Event_value>(-1);

  /**
   * @brief Constructor.
   *
   * @param birth Birth value of the cycle. Default value: @ref inf.
   * @param death Death value of the cycle. Default value: @ref inf.
   * @param dim Dimension of the cycle. Default value: -1.
   */
  Persistence_interval(Event_value birth = inf, Event_value death = inf, Dimension dim = -1)
      : dim(dim), birth(birth), death(death) {}

  Dimension dim;     /**< Dimension of the cycle.*/
  Event_value birth; /**< Birth value of the cycle. */
  Event_value death; /**< Death value of the cycle. */

  /**
   * @brief operator<<
   *
   * @param stream outstream
   * @param interval interval to stream
   */
  inline friend std::ostream& operator<<(std::ostream& stream, const Persistence_interval& interval) {
    stream << "[" << interval.dim << "] ";
    if constexpr (std::numeric_limits<Event_value>::has_infinity) {
      stream << interval.birth << " - " << interval.death;
    } else {
      if (interval.birth == inf)
        stream << "inf";
      else
        stream << interval.birth;
      stream << " - ";
      if (interval.death == inf)
        stream << "inf";
      else
        stream << interval.death;
    }
    return stream;
  }

  /**
   * @brief Specialization of `get` for @ref Gudhi::persistence_matrix::Persistence_interval.
   *
   * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
   * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
   */
  template <std::size_t I>
  constexpr auto& get() & noexcept {
    static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

    if constexpr (I == 0) return birth;
    if constexpr (I == 1) return death;
    if constexpr (I == 2) return dim;
  }

  /**
   * @brief Specialization of `get` for @ref Gudhi::persistence_matrix::Persistence_interval.
   *
   * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
   * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
   */
  template <std::size_t I>
  constexpr const auto& get() const& noexcept {
    static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

    if constexpr (I == 0) return birth;
    if constexpr (I == 1) return death;
    if constexpr (I == 2) return dim;
  }

  /**
   * @brief Specialization of `get` for @ref Gudhi::persistence_matrix::Persistence_interval.
   *
   * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
   * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
   */
  template <std::size_t I>
  constexpr auto&& get() && noexcept {
    static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

    if constexpr (I == 0) return std::move(birth);
    if constexpr (I == 1) return std::move(death);
    if constexpr (I == 2) return std::move(dim);
  }

  /**
   * @brief Specialization of `get` for @ref Gudhi::persistence_matrix::Persistence_interval.
   *
   * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
   * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
   */
  template <std::size_t I>
  constexpr const auto&& get() const&& noexcept {
    static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

    if constexpr (I == 0) return std::move(birth);
    if constexpr (I == 1) return std::move(death);
    if constexpr (I == 2) return std::move(dim);
  }
};

}  // namespace persistence_matrix
}  // namespace Gudhi

namespace std {

/**
 * @ingroup persistence_matrix
 *
 * @brief Partial specialization of `std::tuple_size` for @ref Gudhi::persistence_matrix::Persistence_interval.
 *
 * @tparam Dimension First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam Event_value Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 */
template <typename Dimension, typename Event_value>
struct tuple_size<Gudhi::persistence_matrix::Persistence_interval<Dimension, Event_value> >
    : integral_constant<size_t, 3> {};

/**
 * @ingroup persistence_matrix
 *
 * @brief Partial specialization of `std::tuple_element` for @ref Gudhi::persistence_matrix::Persistence_interval.
 *
 * @tparam I Index of the type to store: 0 for the birth value type, 1 for the death value type and 2 for the
 * dimension value type.
 * @tparam Dimension First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam Event_value Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 */
template <size_t I, typename Dimension, typename Event_value>
struct tuple_element<I, Gudhi::persistence_matrix::Persistence_interval<Dimension, Event_value> > {
  static_assert(I < 3, "Value mismatch at argument 1 in template parameter list. Maximal possible value is 2.");

  using type = typename conditional <I < 2, Event_value, Dimension>::type;
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Partial specialization of `get` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
 * @tparam Dimension First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam Event_value Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @param i Interval from which the value should be returned.
 * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
 */
template <size_t I, typename Dimension, typename Event_value>
constexpr auto& get(Gudhi::persistence_matrix::Persistence_interval<Dimension, Event_value>& i) noexcept {
  return i.template get<I>();
}

/**
 * @ingroup persistence_matrix
 *
 * @brief Partial specialization of `get` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
 * @tparam Dimension First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam Event_value Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @param i Interval from which the value should be returned.
 * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
 */
template <size_t I, typename Dimension, typename Event_value>
constexpr const auto& get(const Gudhi::persistence_matrix::Persistence_interval<Dimension, Event_value>& i) noexcept {
  return i.template get<I>();
}

/**
 * @ingroup persistence_matrix
 *
 * @brief Partial specialization of `get` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
 * @tparam Dimension First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam Event_value Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @param i Interval from which the value should be returned.
 * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
 */
template <size_t I, typename Dimension, typename Event_value>
constexpr auto&& get(Gudhi::persistence_matrix::Persistence_interval<Dimension, Event_value>&& i) noexcept {
  return std::move(i).template get<I>();
}

/**
 * @ingroup persistence_matrix
 *
 * @brief Partial specialization of `get` for @ref Gudhi::persistence_matrix::Persistence_interval.
 * 
 * @tparam I Index of the value to return: 0 for the birth value, 1 for the death value and 2 for the dimension.
 * @tparam Dimension First template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @tparam Event_value Second template parameter of @ref Gudhi::persistence_matrix::Persistence_interval.
 * @param i Interval from which the value should be returned.
 * @return Either the birth value if @p I == 0, the death value if @p I == 1 or the dimension if @p I == 2.
 */
template <size_t I, typename Dimension, typename Event_value>
constexpr const auto&& get(const Gudhi::persistence_matrix::Persistence_interval<Dimension, Event_value>&& i) noexcept {
  return std::move(i).template get<I>();
}

}  // namespace std

#endif  // PM_INTERVAL_INCLUDED
