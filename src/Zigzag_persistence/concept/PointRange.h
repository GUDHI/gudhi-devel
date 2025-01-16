/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ZZ_POINT_RANGE_H_
#define CONCEPT_ZZ_POINT_RANGE_H_

/** @file PointRange.h
 * @brief Contains @ref Gudhi::zigzag_persistence::Point and @ref Gudhi::zigzag_persistence::PointRange concept.
 */

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief Data structure representing a point of fixed dimension. The structure of the point does not matter
 * it-self as long as it corresponds to the input type of the @ref DistanceFunction concept.
 */
class Point{};

/**
 * @brief Range of @ref Point. If used with @ref Oscillating_rips_edge_order_policy::FARTHEST_POINT_ORDERING
 * order policy, it has to be a random access range.
 */
class PointRange {
 public:
  /**
   * @brief Returns begin iterator.
   */
  auto begin();

  /**
   * @brief Returns end iterator.
   */
  auto end();

  /**
   * @brief Returns size of the range.
   */
  std::size_t size();

  /**
   * @brief Necessary only if used with @ref Oscillating_rips_edge_order_policy::FARTHEST_POINT_ORDERING.
   * Returns the element at the given index.
   * 
   * @param index Index of the element to return.
   * @return Point at index @p index.
   */
  Point operator[](std::size_t index);
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // CONCEPT_ZZ_POINT_RANGE_H_
