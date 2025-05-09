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

#include <nanobind/ndarray.h>

using Sequence1D = std::vector<double>;
using Tensor1D = nanobind::ndarray<double, nanobind::ndim<1>>;
using Sequence2D = std::vector<std::vector<double>>;
using Tensor2D = nanobind::ndarray<double, nanobind::ndim<2>>;

inline Sequence2D _get_sequence_from_tensor(const Tensor2D& tensor)
{
  Sequence2D vec(tensor.shape(0), Sequence1D(tensor.shape(1)));
  for (std::size_t i = 0; i < tensor.shape(0); ++i) {
    for (std::size_t j = 0; j < tensor.shape(1); ++j) vec[i][j] = tensor(i, j);
  }
  return vec;
}

inline Sequence1D _get_sequence_from_tensor(const Tensor1D& tensor)
{
  Sequence1D vec(tensor.shape(0));
  for (std::size_t i = 0; i < tensor.shape(0); ++i) {
    vec[i] = tensor(i);
  }
  return vec;
}

using Nearest_landmark_sequence = std::vector<std::vector<std::pair<std::size_t, double>>>;
using Nearest_landmark_tensor = nanobind::ndarray<double, nanobind::shape<-1, -1, 2>>;

inline Nearest_landmark_sequence _get_sequence_from_tensor(const Nearest_landmark_tensor& tensor)
{
  Nearest_landmark_sequence vec(tensor.shape(0), std::vector<std::pair<std::size_t, double>>(tensor.shape(1)));
  for (std::size_t i = 0; i < tensor.shape(0); ++i) {
    for (std::size_t j = 0; j < tensor.shape(1); ++j) {
      vec[i][j].first = tensor(i, j, 0);
      vec[i][j].second = tensor(i, j, 1);
    }
  }
  return vec;
}
