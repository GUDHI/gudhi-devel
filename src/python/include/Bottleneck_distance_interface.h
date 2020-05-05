/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
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

#include <boost/range/iterator_range.hpp>

#include <iostream>
#include <vector>
#include <utility>  // for std::pair

namespace Gudhi {

namespace persistence_diagram {

  // bottleneck_distance function renamed for the python function can be called bottleneck_dstance
  static double bottleneck(void* data1, std::ptrdiff_t n1,
                           void* data2, std::ptrdiff_t n2,
                           double e) {
    auto p1 = static_cast<std::pair<double, double>*>(data1);
    auto p2 = static_cast<std::pair<double, double>*>(data2);
    auto diag1 = boost::make_iterator_range(p1, p1 + n1);
    auto diag2 = boost::make_iterator_range(p2, p2 + n2);
    return bottleneck_distance(diag1, diag2, e);
  }

  static double bottleneck(void* data1, std::ptrdiff_t n1,
                           void* data2, std::ptrdiff_t n2) {
    auto p1 = static_cast<std::pair<double, double>*>(data1);
    auto p2 = static_cast<std::pair<double, double>*>(data2);
    auto diag1 = boost::make_iterator_range(p1, p1 + n1);
    auto diag2 = boost::make_iterator_range(p2, p2 + n2);
    return bottleneck_distance(diag1, diag2);
  }

}  // namespace persistence_diagram

}  // namespace Gudhi


#endif  // INCLUDE_BOTTLENECK_DISTANCE_INTERFACE_H_
