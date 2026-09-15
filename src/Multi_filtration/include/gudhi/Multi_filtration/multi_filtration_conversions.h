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
 * @file multi_filtration_conversions.h
 * @author Hannah Schreiber
 */

#ifndef MF_CONVERSIONS_H_
#define MF_CONVERSIONS_H_

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <limits>
#include <cmath>

#include <boost/safe_numerics/safe_integer.hpp>

#include <gudhi/Multi_filtration/Flat_array_filtration.h>
#include <gudhi/Multi_filtration/Nested_array_filtration.h>
#include <gudhi/Multi_filtration/Degree_bifiltration.h>

namespace Gudhi {

namespace multi_filtration {

/**
 * @ingroup multi_filtration
 *
 * @brief Converts the given @ref Flat_array_filtration into the @ref StoragePolicy given as template argument.
 * It is assumed that the given value is simplified (i.e. minimal and ordered lexicographically). It is also assumed
 * that the given value has the right format for the chosen output type for de default constructor. E.g.,
 * if `OutStoragePolicy` is @ref Degree_bifiltration, then the number of parameters are 2, shift/step are assumed
 * to be 0/1 etc.
 *
 * @tparam OutStoragePolicy New @ref StoragePolicy. Has to be either
 * @ref Gudhi::multi_filtration::Flat_array_filtration,
 * @ref Gudhi::multi_filtration::Nested_array_filtration or
 * @ref Gudhi::multi_filtration::Degree_bifiltration with desired template argument.
 * @tparam Co Only token into account for @ref Gudhi::multi_filtration::Degree_bifiltration. Set to `true` if the
 * filtration value should be interpreted as negative cones instead of positive ones. Default value: false.
 * @tparam T Template parameter of @ref Flat_array_filtration.
 * @param f Filtration value to convert.
 */
template <class OutStoragePolicy, bool Co = false, typename T>
inline OutStoragePolicy as_type(const Flat_array_filtration<T>& f) {
  using U = typename OutStoragePolicy::value_type;

  if constexpr (std::is_same_v<OutStoragePolicy, Flat_array_filtration<U> > ||
                std::is_same_v<OutStoragePolicy, Nested_array_filtration<U> >) {
    const auto& values = f.get_underlying_container();
    return OutStoragePolicy(values.begin(), values.end(), f.num_parameters());
  } else if constexpr (std::is_same_v<OutStoragePolicy, Degree_bifiltration<U> >) {
    if (f.num_parameters() != 2) throw std::invalid_argument("Cannot convert a non-bifiltration to a bifiltration.");

    T maxIndex = 0;
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      maxIndex = maxIndex < f(g, 1) ? f(g, 1) : maxIndex;
    }

    std::vector<U> values(maxIndex + 1, Co ? OutStoragePolicy::template T_m_inf<> : OutStoragePolicy::template T_inf<>);
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      values[f(g, 1)] = f(g, 0);
    }
    return OutStoragePolicy(std::move(values), 2);
  } else {
    throw std::invalid_argument("Given out multi filtration value is not available.");
  }
}

/**
 * @ingroup multi_filtration
 *
 * @brief Converts the given @ref Nested_array_filtration into the @ref StoragePolicy given as template argument.
 * It is assumed that the given value is simplified (i.e. minimal and ordered lexicographically). It is also assumed
 * that the given value has the right format for the chosen output type for de default constructor. E.g.,
 * if `OutStoragePolicy` is @ref Degree_bifiltration, then the number of parameters are 2, shift/step are assumed
 * to be 0/1 etc.
 *
 * @tparam OutStoragePolicy New @ref StoragePolicy. Has to be either
 * @ref Gudhi::multi_filtration::Flat_array_filtration,
 * @ref Gudhi::multi_filtration::Nested_array_filtration or
 * @ref Gudhi::multi_filtration::Degree_bifiltration with desired template argument.
 * @tparam Co Only token into account for @ref Gudhi::multi_filtration::Degree_bifiltration. Set to `true` if the
 * filtration value should be interpreted as negative cones instead of positive ones. Default value: false.
 * @tparam T Template parameter of @ref Nested_array_filtration.
 * @param f Filtration value to convert.
 */
template <class OutStoragePolicy, bool Co = false, typename T>
inline OutStoragePolicy as_type(const Nested_array_filtration<T>& f) {
  using U = typename OutStoragePolicy::value_type;

  if constexpr (std::is_same_v<OutStoragePolicy, Flat_array_filtration<U> >) {
    std::vector<U> values(f.num_entries());
    std::size_t i = 0;
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      for (std::size_t p = 0; p < f.num_parameters(); ++p) {
        values[i] = f(g, p);
        ++i;
      }
    }
    return OutStoragePolicy(std::move(values), f.num_parameters());
  } else if constexpr (std::is_same_v<OutStoragePolicy, Nested_array_filtration<U> >) {
    std::vector<std::vector<U> > out(f.num_generators());
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      out[g] = std::vector<U>(f.begin(g), f.end(g));
    }
    return Nested_array_filtration<U>(std::move(out), f.num_parameters());
  } else if constexpr (std::is_same_v<OutStoragePolicy, Degree_bifiltration<U> >) {
    if (f.num_parameters() != 2) throw std::invalid_argument("Cannot convert a non-bifiltration to a bifiltration.");

    T maxIndex = 0;
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      maxIndex = maxIndex < f(g, 1) ? f(g, 1) : maxIndex;
    }

    std::vector<U> values(maxIndex + 1, Co ? OutStoragePolicy::template T_m_inf<> : OutStoragePolicy::template T_inf<>);
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      values[f(g, 1)] = f(g, 0);
    }
    return OutStoragePolicy(std::move(values), 2);
  } else {
    throw std::invalid_argument("Given out multi filtration value is not available.");
  }
}

/**
 * @ingroup multi_filtration
 *
 * @brief Converts the given @ref Degree_bifiltration into the @ref StoragePolicy given as template argument.
 *
 * @tparam OutStoragePolicy New @ref StoragePolicy. Has to be either
 * @ref Gudhi::multi_filtration::Flat_array_filtration,
 * @ref Gudhi::multi_filtration::Nested_array_filtration or
 * @ref Gudhi::multi_filtration::Degree_bifiltration with desired template argument.
 * @tparam Co Only token into account for @ref Gudhi::multi_filtration::Flat_array_filtration and
 * @ref Gudhi::multi_filtration::Nested_array_filtration. Set to `true` if the filtration value should be
 * interpreted as negative cones instead of positive ones. Default value: false.
 * @tparam T Template parameter of @ref Degree_bifiltration.
 * @param f Filtration value to convert.
 */
template <class OutStoragePolicy, bool Co = false, typename T>
inline OutStoragePolicy as_type(const Degree_bifiltration<T>& f) {
  using U = typename OutStoragePolicy::value_type;

  if constexpr (std::is_same_v<OutStoragePolicy, Degree_bifiltration<U> >) {
    std::vector<U> out(f.get_underlying_container().begin(), f.get_underlying_container().end());
    Degree_bifiltration<U> res(std::move(out), 2);
    res.set_mapping(f.get_shift(), f.get_step());
    return res;
  } else {
    auto gen_index = [&f](std::size_t i) {
      if constexpr (Co) {
        return f.num_generators() - 1 - i;
      } else {
        return i;
      }
    };

    auto strictly_dominates = [](T a, T b) {
      if constexpr (Co) {
        return a < b;
      } else {
        return a > b;
      }
    };

    if (f.num_generators() == 0) return OutStoragePolicy();

    std::vector<std::size_t> order(f.num_generators());
    std::iota(order.begin(), order.end(), 0);
    // lexicographical order
    std::sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j) {
      if (f(i, 0) == f(j, 0)) return f(i, 1) < f(j, 1);  // f(i, 1) and f(j, 1) cannot be equal for i != j
      return f(i, 0) < f(j, 0);
    });

    if constexpr (std::is_same_v<OutStoragePolicy, Flat_array_filtration<U> >) {
      std::vector<U> values;
      values.reserve(f.num_generators() * 2);
      std::size_t g = order[gen_index(0)];
      T threshold = f(g, 1);
      values.push_back(f(g, 0));
      values.push_back(threshold);
      for (std::size_t i = 1; i < f.num_generators(); ++i) {
        g = order[gen_index(i)];
        if (strictly_dominates(threshold, f(g, 1))) {
          threshold = f(g, 1);
          values.push_back(f(g, 0));
          values.push_back(threshold);
        }
      }
      if constexpr (Co) {
        // lexicographical order
        const std::size_t max_idx = values.size() - 1;
        for (std::size_t i = 0; i < values.size() / 2; i += 2) {
          std::swap(values[i], values[max_idx - 1 - i]);
          std::swap(values[i + 1], values[max_idx - i]);
        }
      }

      return OutStoragePolicy(std::move(values), 2);
    } else if constexpr (std::is_same_v<OutStoragePolicy, Nested_array_filtration<U> >) {
      std::vector<std::vector<U> > values;
      values.reserve(f.num_generators());
      std::size_t g = order[gen_index(0)];
      U threshold = f(g, 1);
      values.emplace_back(std::vector<U>{static_cast<U>(f(g, 0)), threshold});
      for (std::size_t i = 1; i < f.num_generators(); ++i) {
        g = order[gen_index(i)];
        if (strictly_dominates(threshold, f(g, 1))) {
          threshold = f(g, 1);
          std::vector<U> v = {static_cast<U>(f(g, 0)), threshold};
          values.emplace_back(std::move(v));
        }
      }
      if constexpr (Co) {
        // lexicographical order
        std::reverse(values.begin(), values.end());
      }

      return OutStoragePolicy(std::move(values), 2);
    } else {
      throw std::invalid_argument("Given out multi filtration value is not available.");
    }
  }
}

namespace detail {

/**
 * @private
 * Tries as good as it can to be precise, but there is probably room for improvements...
 */
template <typename T, typename U>
bool _safe_equal(T t, U u) {
  static_assert(std::is_arithmetic_v<T> && std::is_arithmetic_v<U>, "_safe_equal requires arithmetic types");

  if constexpr (std::is_same_v<T, U>) {
    return t == u;
  } else if constexpr (std::is_integral_v<T> && std::is_integral_v<U>) {
    // integer x integer
    // in C++20, can be replaced with std::cmp_equal
    return boost::safe_numerics::safe<T>(t) == boost::safe_numerics::safe<U>(u);
  } else if constexpr (std::is_floating_point_v<T> && std::is_floating_point_v<U>) {
    // floating x floating
    using Wider = std::conditional_t<(sizeof(T) >= sizeof(U)), T, U>;
    return static_cast<Wider>(t) == static_cast<Wider>(u);
  } else {
    if constexpr (std::is_floating_point_v<T>) {
      // floating x integer
      // if the floating type is nan or has non-zero decimal values, it cannot be equal to an integer.
      if (std::isnan(t) || std::trunc(t) != t) return false;
      // the floating type does not have any decimal, i.e. is comparable with an integer
      if (t < static_cast<T>(std::numeric_limits<U>::min()) || t > static_cast<T>(std::numeric_limits<U>::max())) {
        return false;
      }
      return static_cast<U>(t) == u;
    } else {
      // integer x floating
      if (std::isnan(u) || std::trunc(u) != u) return false;
      if (u < static_cast<U>(std::numeric_limits<T>::min()) || u > static_cast<U>(std::numeric_limits<T>::max())) {
        return false;
      }
      return t == static_cast<T>(u);
    }
  }
}

}  // namespace detail

/**
 * @brief Compares for equality two filtration values potentially not of the same type.
 *
 * @tparam MultiFiltrationValue1 A filtration value type with methods: num_parameters(), num_generators(), is_nan()
 * and operator(g, p).
 * @tparam MultiFiltrationValue2 A filtration value type with methods: num_parameters(), num_generators(), is_nan()
 * and operator(g, p).
 */
template <class MultiFiltrationValue1, class MultiFiltrationValue2>
inline bool are_equal_filtration_values(const MultiFiltrationValue1& f1, const MultiFiltrationValue2& f2) {
  if constexpr (std::is_same_v<MultiFiltrationValue1, MultiFiltrationValue2>) {
    if (&f1 == &f2) return true;
  }
  if (f1.is_nan() || f2.is_nan()) return false;
  if (f1.num_parameters() != f2.num_parameters()) return false;
  if (f1.num_generators() != f2.num_generators()) return false;

  for (std::size_t p = 0; p < f1.num_parameters(); ++p) {
    for (std::size_t g = 0; g < f1.num_generators(); ++g) {
      if (!detail::_safe_equal(f1(g, p), f2(g, p))) return false;
    }
  }
  return true;
}

}  // namespace multi_filtration

}  // namespace Gudhi

#endif  // MF_CONVERSIONS_H_
