/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 INRIA
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
