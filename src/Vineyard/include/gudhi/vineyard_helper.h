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
#include <stdexcept>
#include <vector>

#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace vineyard {

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

  Index numberOfSimplices = boundaries.size();

  for (auto sh : complex.complex_simplex_range()) {
    complex.assign_key(sh, numberOfSimplices);
    ++numberOfSimplices;
  }

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

}  // namespace vineyard
}  // namespace Gudhi

#endif  // GUDHI_VINEYARD_HELPER_H_
