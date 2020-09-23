/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATH_H_
#define MATH_H_

#include <boost/math/constants/constants.hpp>

namespace Gudhi {

// In wait of C++20 std::numbers::pi with #include <numbers>
static constexpr double PI = boost::math::constants::pi<double>();

}  // namespace Gudhi

#endif   // MATH_H_
