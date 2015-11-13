/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
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
#ifndef DEBUG_UTILS_H_
#define DEBUG_UTILS_H_

#include <iostream>

#ifdef NDEBUG
  // GUDHI_NDEBUG is the Gudhi official flag for debug mode.
  #define GUDHI_NDEBUG
#endif

#define PRINT(a) std::cerr << #a << ": " << (a) << " (DISP)" << std::endl

#ifdef GUDHI_NDEBUG
  #define DBG(a) std::cout << "DBG: " << (a) << std::endl
  #define DBGMSG(a, b) std::cout << "DBG: " << a << b << std::endl
  #define DBGVALUE(a) std::cout << "DBG: " <<  #a << ": " << a << std::endl
  #define DBGCONT(a) std::cout << "DBG: container " << #a << " -> "; for (auto x : a) std::cout << x << ","; std::cout << std::endl
#else
  #define DBG(a) (void) 0
  #define DBGMSG(a, b) (void) 0
  #define DBGVALUE(a) (void) 0
  #define DBGCONT(a) (void) 0
#endif

// GUDHI_CHECK throw an exception on condition in debug mode, but does nothing in release mode
#ifdef GUDHI_NDEBUG
  #define GUDHI_CHECK(cond, excpt) if (cond) throw excpt
#else
  #define GUDHI_CHECK(cond, excpt) (void) 0
#endif

#endif  // DEBUG_UTILS_H_
