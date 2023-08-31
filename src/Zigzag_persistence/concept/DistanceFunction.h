/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ZZ_DISTANCE_FUNCTION_H_
#define CONCEPT_ZZ_DISTANCE_FUNCTION_H_

/** @file DistanceFunction.h
 * @brief Contains @ref Gudhi::zigzag_persistence::DistanceFunction concept.
 */

#include "PointRange.h"

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief Distance function taking two points as input and returning the distance between them.
 */
using DistanceFunction = double (*)(const Point&, const Point&);

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // CONCEPT_ZZ_DISTANCE_FUNCTION_H_
