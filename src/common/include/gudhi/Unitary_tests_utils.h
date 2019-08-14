/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2017 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
#ifndef UNITARY_TESTS_UTILS_H_
#define UNITARY_TESTS_UTILS_H_

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <limits>  // for std::numeric_limits<>

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
