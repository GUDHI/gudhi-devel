/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

// #include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

using Sequence = std::vector<std::vector<double>>;
using Tensor = nanobind::ndarray<double, nanobind::ndim<2>>;

inline Sequence _get_sequence_from_tensor(const Tensor& tensor)
{
  Sequence vec(tensor.shape(0), std::vector<double>(tensor.shape(1)));
  for (std::size_t i = 0; i < tensor.shape(0); ++i) {
    for (std::size_t j = 0; j < tensor.shape(1); ++j) vec[i][j] = tensor(i, j);
  }
  return vec;
}
