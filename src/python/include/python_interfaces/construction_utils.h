/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_INCLUDE_CONSTRUCTION_UTILS_PYTHON_H_
#define GUDHI_INCLUDE_CONSTRUCTION_UTILS_PYTHON_H_

#include <cstddef>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "numpy_utils.h"

// TODO: add this namespace to other helper files? They are only used in private anyway.
namespace Gudhi {
namespace python {

template <typename T>
static nanobind::object _to_nanobind_object(T &&val) {
  using U = std::decay_t<T>;
  // to avoid non trivial nanobind::cast when not necessary
  if constexpr (std::is_base_of_v<nanobind::object, U>) {
    return std::forward<T>(val);
  } else {
    return nanobind::cast(std::forward<T>(val));
  }
}

// more efficient than a list.append if n is known in advance
// (and the object does not have to be resized later obviously)
template <class F>
nanobind::tuple _build_tuple(std::size_t n, F &&construct_row) {
  PyObject *out = PyTuple_New(static_cast<Py_ssize_t>(n));
  if (!out) throw nanobind::python_error();

  for (std::size_t p = 0; p < n; ++p) {
    nanobind::object row = _to_nanobind_object(std::forward<F>(construct_row)(p));
    PyTuple_SET_ITEM(out, p, row.release().ptr());
  }

  return nanobind::steal<nanobind::tuple>(out);
}

template <typename Integer, typename Arithmetic,
          class = std::enable_if_t<std::is_integral_v<Integer> && std::is_arithmetic_v<Arithmetic>>>
bool _is_and_fits_in_int(Arithmetic x) {
  using Out_limits = std::numeric_limits<Integer>;

  if constexpr (std::is_floating_point_v<Arithmetic>) {
    return std::isfinite(x) && std::trunc(x) == x && x >= static_cast<Arithmetic>(Out_limits::min()) &&
           x <= static_cast<Arithmetic>(Out_limits::max());
  } else {
    if constexpr (std::is_signed_v<Arithmetic> == std::is_signed_v<Integer>) {
      return x >= Out_limits::min() && x <= Out_limits::max();
    } else if constexpr (std::is_signed_v<Arithmetic>) {
      if (x < 0) return false;
      return static_cast<std::make_unsigned_t<Arithmetic>>(x) <= Out_limits::max();
    } else {
      using Wider_unsigned =
          std::make_unsigned_t<std::conditional_t<(sizeof(Arithmetic) > sizeof(Integer)), Arithmetic, Integer>>;
      return static_cast<Wider_unsigned>(x) <= static_cast<Wider_unsigned>(Out_limits::max());
    }
  }
}

template <typename Integer, typename Arithmetic,
          class = std::enable_if_t<std::is_integral_v<Integer> && std::is_arithmetic_v<Arithmetic>>>
Integer _cast_to_int(Arithmetic x,
                     const std::string &errorMsg = "Given value does not represent or fits into an integer.") {
  if (!_is_and_fits_in_int<Integer>(x)) throw std::runtime_error(errorMsg.c_str());
  return static_cast<Integer>(x);
}

template <typename T>
struct is_ndarray : std::false_type {};

template <typename T, typename... Args>
struct is_ndarray<nanobind::ndarray<T, Args...>> : std::true_type {};

template <typename T>
struct contains_ndarray : is_ndarray<T> {};

// works only for std::vector, but if other containers become necessary, they can be added the same way
template <typename T>
struct contains_ndarray<std::vector<T>> : contains_ndarray<T> {};

template <typename T>
inline constexpr bool contains_ndarray_v = contains_ndarray<T>::value;

template <typename Variant, typename T, typename... Rest>
static bool _try_fill(nanobind::iterable data, Variant &out) {
  // Avoids atemps to copy the data if it contains an ndarray
  constexpr bool convert = !contains_ndarray_v<T>;

  if (T val; nanobind::try_cast<T>(data, val, convert)) {
    out = std::move(val);
    return true;
  }

  if constexpr (sizeof...(Rest) > 0)
    return _try_fill<Variant, Rest...>(data, out);
  else
    return false;
}

// Carefull, order in Ts matters! Use types which do not need copies first, such that they are tested first!
template <typename... Ts>
static auto _convert_iterable_to_cpp_type(
    nanobind::iterable data, const std::string &errorMessage = "Iterable does not match any of the accepted types.") {
  static_assert(sizeof...(Ts) > 0, "Needs at least one candidate type.");
  using Variant = std::variant<Ts...>;

  Variant result{};  // the first alternative has to be default-constructable
  if (!_try_fill<Variant, Ts...>(data, result)) {
    throw nanobind::type_error(errorMessage.c_str());
  }
  return result;
}

template <typename T>
struct wrap_ndarray {
  using type = T;
};

template <typename T>
struct wrap_ndarray<nanobind::ndarray<const T, nanobind::ndim<1>, nanobind::any_contig>> {
  using type = Numpy_span<T>;
};

template <typename T>
struct wrap_ndarray<nanobind::ndarray<const T, nanobind::ndim<2>>> {
  using type = Numpy_2d_span<T>;
};

// works only for std::vector, but if other containers become necessary, they can be added the same way
template <typename T>
struct wrap_ndarray<std::vector<T>> {
  using type = std::vector<typename wrap_ndarray<T>::type>;
};

template <typename T>
using wrap_ndarray_t = typename wrap_ndarray<T>::type;

template <typename T>
wrap_ndarray_t<std::decay_t<T>> _wrap_value(T &&val) {
  using U = std::decay_t<T>;

  if constexpr (!contains_ndarray_v<U>) {
    // nothing to do
    return std::forward<T>(val);
  } else if constexpr (is_ndarray<U>::value) {
    return wrap_ndarray_t<U>(std::forward<T>(val));
  } else {
    // U is e.g. std::vector<ndarray<...>> or std::vector<std::vector<ndarray<...>>>
    wrap_ndarray_t<U> out;
    out.reserve(val.size());
    for (auto &elem : val) out.push_back(_wrap_value(std::move(elem)));
    return out;
  }
}

template <typename Variant, typename T, typename... Rest>
static bool _try_fill_and_wrap_ndarrays(nanobind::iterable data, std::optional<Variant> &out) {
  constexpr bool convert = !contains_ndarray_v<T>;

  if (T val; nanobind::try_cast<T>(data, val, convert)) {
    out.emplace(std::in_place_type<wrap_ndarray_t<std::decay_t<T>>>, _wrap_value(std::move(val)));
    return true;
  }

  if constexpr (sizeof...(Rest) > 0)
    return _try_fill_and_wrap_ndarrays<Variant, Rest...>(data, out);
  else
    return false;
}

// Carefull, order in Ts matters! Use types which do not need copies first, such that they are tested first!
template <typename... Ts>
static auto _convert_iterable_to_cpp_type_and_wrap_ndarrays(
    nanobind::iterable data, const std::string &errorMessage = "Iterable does not match any of the accepted types.") {
  static_assert(sizeof...(Ts) > 0, "Needs at least one candidate type.");
  using Variant = std::variant<wrap_ndarray_t<Ts>...>;

  std::optional<Variant> result;  // optional to avoid obligation of the first alternative to be default-constructable
  if (!_try_fill_and_wrap_ndarrays<Variant, Ts...>(data, result)) {
    throw nanobind::type_error(errorMessage.c_str());
  }
  return std::move(*result);
}

}  // namespace python
}  // namespace Gudhi

#endif  // GUDHI_INCLUDE_CONSTRUCTION_UTILS_PYTHON_H_
