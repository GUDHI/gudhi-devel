/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_BOTTLENECK_DISTANCE_INTERFACE_H_
#define INCLUDE_BOTTLENECK_DISTANCE_INTERFACE_H_

#include <gudhi/Bottleneck.h>

#include <iostream>
#include <vector>
#include <utility>  // for std::pair

namespace Gudhi {

namespace persistence_diagram {

  // bottleneck_distance function renamed for the python function can be called bottleneck_dstance
  double bottleneck(const std::vector<std::pair<double, double>>& diag1,
                    const std::vector<std::pair<double, double>>& diag2,
                    double e) {
    return bottleneck_distance(diag1, diag2, e);
  }

  double bottleneck(const std::vector<std::pair<double, double>>& diag1,
                    const std::vector<std::pair<double, double>>& diag2) {
    return bottleneck_distance(diag1, diag2);
  }

}  // namespace persistence_diagram

}  // namespace Gudhi


#endif  // INCLUDE_BOTTLENECK_DISTANCE_INTERFACE_H_
