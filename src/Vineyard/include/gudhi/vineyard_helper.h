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
 * @ingroup vineyard
 * 
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

/**
 * @ingroup vineyard
 * 
 * @brief Builds a boundary matrix from the given matrix.
 * 
 * @tparam FilteredComplex Complex type with following methods (see @ref FilteredComplex or @ref Gudhi::Simplex_tree
 * for method description): complex_simplex_range, skeleton_simplex_range, key, assign_key, dimension, filtration
 * and boundary_simplex_range. And the following types: Filtration_value and Simplex_key.
 * @tparam Filtration_value Filtration value type for the boundary matrix. `FilteredComplex::Filtration_value` has
 * to be convertible into it. Default value: FilteredComplex::Filtration_value.
 * @tparam Index Index type for the boundary matrix. Has to be an integer type and big enough to count all cells in
 * the complex. Default value: FilteredComplex::Simplex_key.
 * @tparam Dimension Dimension type. Has to be an integer type.
 * @param[in] complex Complex to convert.
 * @param[out] boundaries Container for the boundaries. If not empty, the elements are not erased and the new
 * boundaries are added at the end.
 * @param[out] dimensions Container for the dimensions. If not empty, the elements are not erased and the new
 * dimensions are added at the end. Has to have the same size than @p boundaries as the indices should correspond.
 * @param[out] filtrationValues Container for the filtration values. If not empty, the elements are not erased and
 * the new values are added at the end. Has to have the same size than @p boundaries as the indices should correspond.
 */
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

/**
 * @ingroup vineyard
 * 
 * @brief Same method than @ref build_boundary_matrix_from_complex(FilteredComplex&, std::vector<std::vector<Index> >&,
 * std::vector<Dimension>&, std::vector<Filtration_value>&), but only the filtration values are computed. The order of
 * the filtration values are the same than if computed with the other version. Usefull if one already called the other
 * version and then only modified the filtration values of the complex and therefore don't need to recompute the
 * boundaries.
 *
 * @tparam FilteredComplex Complex type with following methods (see @ref FilteredComplex or @ref Gudhi::Simplex_tree
 * for method description): complex_simplex_range, skeleton_simplex_range, key, assign_key, dimension, filtration
 * and boundary_simplex_range. And the following types: Filtration_value and Simplex_key.
 * @tparam Filtration_value Filtration value type for the boundary matrix. `FilteredComplex::Filtration_value` has
 * to be convertible into it. Default value: FilteredComplex::Filtration_value.
 * @param[in] complex Complex to convert.
 * @param[out] filtrationValues Container for the filtration values. If not empty, the elements are not erased and
 * the new values are added at the end.
 */
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
