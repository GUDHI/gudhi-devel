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

#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>
#include <gudhi/Degree_rips_bifiltration.h>

namespace Gudhi {

namespace multi_filtration {

/**
 * @brief Converts the given multi filtration value into the type given as template argument. It is assumed that the
 * given value is simplified (i.e. minimal and ordered lexicographically). If the new type is
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration it is additionally assumed that the given value is compatible
 *  with the type, that is,the number of parameters is 2 and the first parameter is an index (positive and convertible
 * to an integer without loss).
 *
 * @tparam Out_multi_filtration New filtration value type. Has to be either
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration with desired template arguments.
 * @tparam T First template argument of the initial filtration value type.
 * @tparam Co Second template argument of the initial filtration value type.
 * @tparam Ensure1Criticality Third template argument of the initial filtration value type.
 * @param f Filtration value to convert.
 */
template <class Out_multi_filtration, typename T, bool Co, bool Ensure1Criticality>
Out_multi_filtration as_type(const Multi_parameter_filtration<T, Co, Ensure1Criticality>& f)
{
  using U = typename Out_multi_filtration::value_type;
  constexpr bool co = Out_multi_filtration::has_negative_cones();
  constexpr bool one_crit = Out_multi_filtration::ensures_1_criticality();

  if constexpr (std::is_same_v<Out_multi_filtration, Multi_parameter_filtration<U, co, one_crit> >) {
    return f.template as_type<U, co, one_crit>();
  } else if constexpr (std::is_same_v<Out_multi_filtration, Dynamic_multi_parameter_filtration<U, co, one_crit> >) {
    return Out_multi_filtration(f.begin(), f.end(), f.num_parameters());
  } else if constexpr (std::is_same_v<Out_multi_filtration, Degree_rips_bifiltration<U, co, one_crit> >) {
    if (f.num_parameters() != 2) throw std::invalid_argument("Cannot convert a non-bifiltration to a bifiltration.");
    U inf = co ? Out_multi_filtration::T_m_inf : Out_multi_filtration::T_inf;

    // f generators assumed to be in lexicographical order
    std::vector<U> values(f(f.num_generators() - 1, 0) + 1, inf);
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      values[f(g, 0)] = f(g, 1);
    }
    return Out_multi_filtration(std::move(values), 2);
  } else {
    throw std::invalid_argument("Given out multi filtration value is not available.");
  }
}

/**
 * @brief Converts the given multi filtration value into the type given as template argument. It is assumed that the
 * given value is simplified (i.e. minimal and ordered lexicographically). If the new type is
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration it is additionally assumed that the given value is compatible
 *  with the type, that is,the number of parameters is 2 and the first parameter is an index (positive and convertible
 * to an integer without loss).
 *
 * @tparam Out_multi_filtration New filtration value type. Has to be either
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration with desired template arguments.
 * @tparam T First template argument of the initial filtration value type.
 * @tparam Co Second template argument of the initial filtration value type.
 * @tparam Ensure1Criticality Third template argument of the initial filtration value type.
 * @param f Filtration value to convert.
 */
template <class Out_multi_filtration, typename T, bool Co, bool Ensure1Criticality>
Out_multi_filtration as_type(const Dynamic_multi_parameter_filtration<T, Co, Ensure1Criticality>& f)
{
  using U = typename Out_multi_filtration::value_type;
  constexpr bool co = Out_multi_filtration::has_negative_cones();
  constexpr bool one_crit = Out_multi_filtration::ensures_1_criticality();

  if constexpr (std::is_same_v<Out_multi_filtration, Multi_parameter_filtration<U, co, one_crit> >) {
    std::vector<U> values(f.num_entries());
    std::size_t i = 0;
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      for (std::size_t p = 0; p < f.num_parameters(); ++p) {
        values[i] = f(g, p);
        ++i;
      }
    }
    return Out_multi_filtration(std::move(values), f.num_parameters());
  } else if constexpr (std::is_same_v<Out_multi_filtration, Dynamic_multi_parameter_filtration<U, co, one_crit> >) {
    return f.template as_type<U, co, one_crit>();
  } else if constexpr (std::is_same_v<Out_multi_filtration, Degree_rips_bifiltration<U, co, one_crit> >) {
    if (f.num_parameters() != 2) throw std::invalid_argument("Cannot convert a non-bifiltration to a bifiltration.");
    U inf = co ? Out_multi_filtration::T_m_inf : Out_multi_filtration::T_inf;

    // f generators assumed to be in lexicographical order
    std::vector<U> values(f(f.num_generators() - 1, 0) + 1, inf);
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      values[f(g, 0)] = f(g, 1);
    }
    return Out_multi_filtration(std::move(values), 2);
  } else {
    throw std::invalid_argument("Given out multi filtration value is not available.");
  }
}

/**
 * @brief Converts the given multi filtration value into the type given as template argument.
 *
 * @tparam Out_multi_filtration New filtration value type. Has to be either
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration with desired template arguments.
 * @tparam T First template argument of the initial filtration value type.
 * @tparam Co Second template argument of the initial filtration value type.
 * @tparam Ensure1Criticality Third template argument of the initial filtration value type.
 * @param f Filtration value to convert.
 */
template <class Out_multi_filtration, typename T, bool Co, bool Ensure1Criticality>
Out_multi_filtration as_type(const Degree_rips_bifiltration<T, Co, Ensure1Criticality>& f)
{
  using U = typename Out_multi_filtration::value_type;
  constexpr bool co = Out_multi_filtration::has_negative_cones();
  constexpr bool one_crit = Out_multi_filtration::ensures_1_criticality();

  if constexpr (std::is_same_v<Out_multi_filtration, Degree_rips_bifiltration<U, co, one_crit> >) {
    return f.template as_type<U, co, one_crit>();
  } else {
    auto gen_index = [&f, co](std::size_t i) {
      if constexpr (co) {
        return f.num_generators() - 1 - i;
      } else {
        return i;
      }
    };

    auto strictly_dominates = [co](T a, T b) {
      if constexpr (co) {
        return a < b;
      } else {
        return a > b;
      }
    };

    if (f.size() == 0) return Out_multi_filtration(0);

    if constexpr (std::is_same_v<Out_multi_filtration, Multi_parameter_filtration<U, co, one_crit> >) {
      std::vector<U> values;
      values.reserve(f.num_generators() * 2);
      T threshold = f(gen_index(0), 1);
      values.push_back(gen_index(0));
      values.push_back(threshold);
      for (std::size_t i = 1; i < f.num_generators(); ++i) {
        std::size_t g = gen_index(i);
        if (strictly_dominates(threshold, f(g, 1))) {
          threshold = f(g, 1);
          values.push_back(g);
          values.push_back(threshold);
        }
      }
      if constexpr (co) {
        // lexicographical order
        const std::size_t max_idx = values.size() - 1;
        for (std::size_t i = 0; i < values.size() / 2; i += 2) {
          std::swap(values[i], values[max_idx - 1 - i]);
          std::swap(values[i + 1], values[max_idx - i]);
        }
      }

      return Out_multi_filtration(std::move(values), 2);
    } else if constexpr (std::is_same_v<Out_multi_filtration, Dynamic_multi_parameter_filtration<U, co, one_crit> >) {
      std::vector<Multi_parameter_generator<U> > values;
      values.reserve(f.num_generators());
      T threshold = f(gen_index(0), 1);
      values.emplace_back(std::vector<T>{static_cast<T>(gen_index(0)), threshold});
      for (std::size_t i = 1; i < f.num_generators(); ++i) {
        std::size_t g = gen_index(i);
        if (strictly_dominates(threshold, f(g, 1))) {
          threshold = f(g, 1);
          std::vector<T> v = {static_cast<T>(g), threshold};
          values.emplace_back(std::move(v));
        }
      }
      if constexpr (co) {
        // lexicographical order
        std::reverse(values.begin(), values.end());
      }

      return Out_multi_filtration(std::move(values), 2);
    } else {
      throw std::invalid_argument("Given out multi filtration value is not available.");
    }
  }
}

}  // namespace multi_filtration

}  // namespace Gudhi

#endif  // MF_CONVERSIONS_H_
