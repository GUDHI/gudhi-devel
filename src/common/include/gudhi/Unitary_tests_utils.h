/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2017 Inria
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
#ifndef UNITARY_TESTS_UTILS_H_
#define UNITARY_TESTS_UTILS_H_

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <limits>  // for std::numeric_limits<>
#include <cmath>

template<typename FloatingType >
void GUDHI_TEST_FLOAT_EQUALITY_CHECK(FloatingType a, FloatingType b,
                                     FloatingType epsilon = std::numeric_limits<FloatingType>::epsilon()) {
#ifdef DEBUG_TRACES
  std::cout << "GUDHI_TEST_FLOAT_EQUALITY_CHECK - " << a << " versus " << b
            << " | diff = " << std::fabs(a - b) << " - epsilon = " << epsilon << std::endl;
#endif
  BOOST_CHECK(std::fabs(a - b) <= epsilon);
}

#endif  // UNITARY_TESTS_UTILS_H_
