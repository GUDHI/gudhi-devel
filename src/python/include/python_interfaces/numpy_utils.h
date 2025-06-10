/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_NUMPY_UTILS_PYTHON_H_
#define INCLUDE_NUMPY_UTILS_PYTHON_H_

#include <cstddef>
#include <vector>

#include <nanobind/ndarray.h>

template <class T, typename... Shape>
inline auto _wrap_as_numpy_array(std::vector<T> *tensor, Shape... shapes)
{
  return nanobind::ndarray<nanobind::numpy, T>(
      tensor->data(), {static_cast<std::size_t>(shapes)...}, nanobind::capsule(tensor, [](void *p) noexcept {
        delete reinterpret_cast<std::vector<T> *>(p);
      }));
}

// tensor has to be declared with 'new []'
template <class T, typename... Shape>
inline auto _wrap_as_numpy_array(T *tensor, Shape... shapes)
{
  return nanobind::ndarray<nanobind::numpy, T>(
      tensor, {static_cast<std::size_t>(shapes)...}, nanobind::capsule(tensor, [](void *p) noexcept {
        delete[] reinterpret_cast<T *>(p);
      }));
}

#endif  // INCLUDE_NUMPY_UTILS_PYTHON_H_
