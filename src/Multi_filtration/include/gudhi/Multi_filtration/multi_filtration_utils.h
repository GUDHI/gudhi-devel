/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MF_UTILS_H_
#define MF_UTILS_H_

#include <cstddef>
#include <type_traits>

namespace Gudhi {

namespace multi_filtration {

template <typename T>
class RangeTraits
{
 private:
  static auto check_begin(...) -> std::false_type;
  template <typename U>
  static auto check_begin(U x) -> decltype(x.begin(), std::true_type{});

  static auto check_dynamic_filtration(...) -> std::false_type;
  template <typename U>
  static auto check_dynamic_filtration(U x) -> decltype(x.operator[](std::size_t{}), std::true_type{});

  static auto check_filtration(...) -> std::false_type;
  template <typename U>
  static auto check_filtration(U x) -> decltype(x.ensures_1_criticality(), std::true_type{});

 public:
  static constexpr bool has_begin = decltype(check_begin(std::declval<T>()))::value;
  static constexpr bool is_dynamic_multi_filtration = decltype(check_filtration(std::declval<T>()))::value &&
                                                      decltype(check_dynamic_filtration(std::declval<T>()))::value;
};

}  // namespace multi_filtration

}  // namespace Gudhi

#endif  // MF_UTILS_H_
