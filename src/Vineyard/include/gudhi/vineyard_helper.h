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
 * @file vineyard_helper.h
 * @author Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::vineyard::build_boundary_matrix_from_complex method.
 */

#ifndef GUDHI_VINEYARD_HELPER_H_
#define GUDHI_VINEYARD_HELPER_H_

#include <algorithm>
#include <vector>

#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace vineyard {

/**
 * @private
 */
template <class FilteredComplex, typename Index>
inline Index assign_keys_(FilteredComplex& complex, Index start)
{
  Index numberOfSimplices = start;
  // Vertex ID should correspond to position in the original point cloud if there was one.
  for (auto sh : complex.skeleton_simplex_range(0)) {
    complex.assign_key(sh, numberOfSimplices);
    ++numberOfSimplices;
  }
  for (auto sh : complex.complex_simplex_range()) {
    if (complex.dimension(sh) != 0) {
      complex.assign_key(sh, numberOfSimplices);
      ++numberOfSimplices;
    }
  }
  return numberOfSimplices;
}

template <class FilteredComplex,
          typename Filtration_value = typename FilteredComplex::Filtration_value,
          typename Index = typename FilteredComplex::Simplex_key,
          typename Dimension = int>
inline void build_boundary_matrix_from_complex(FilteredComplex& complex,
                                               std::vector<std::vector<Index> >& boundaries,
                                               std::vector<Dimension>& dimensions,
                                               std::vector<Filtration_value>& filtrationValues)
{
  GUDHI_CHECK(boundaries.size() == dimensions.size() && boundaries.size() == filtrationValues.size(),
              std::invalid_argument("Output containers do not start with the same size"));

  Index numberOfSimplices = assign_keys_(complex, boundaries.size());

  boundaries.resize(numberOfSimplices);
  dimensions.resize(numberOfSimplices);
  filtrationValues.resize(numberOfSimplices);

  for (auto sh : complex.complex_simplex_range()) {
    auto index = complex.key(sh);
    dimensions[index] = complex.dimension(sh);
    filtrationValues[index] = complex.filtration(sh);
    std::vector<Index> boundary;
    for (auto b : complex.boundary_simplex_range(sh)) {
      boundary.push_back(complex.key(b));
    }
    std::sort(boundary.begin(), boundary.end());
    boundaries[index] = std::move(boundary);
  }
}

// same name to emphasize that the filtration values computed are in exactly the same order than if computed 
// with the other version. Usefull if you already called the other version and then only modified the filtration
// values of the complex and you don't need to recompute the boundaries.
template <class FilteredComplex, typename Filtration_value = typename FilteredComplex::Filtration_value>
inline void build_boundary_matrix_from_complex(FilteredComplex& complex,
                                               std::vector<Filtration_value>& filtrationValues)
{
  auto numberOfSimplices = assign_keys_(complex, filtrationValues.size());

  filtrationValues.resize(numberOfSimplices);

  for (auto sh : complex.complex_simplex_range()) {
    filtrationValues[complex.key(sh)] = complex.filtration(sh);
  }
}

}  // namespace vineyard
}  // namespace Gudhi

#endif  // GUDHI_VINEYARD_HELPER_H_
