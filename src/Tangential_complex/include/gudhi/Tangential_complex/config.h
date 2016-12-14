/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
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

#ifndef TANGENTIAL_COMPLEX_CONFIG_H_
#define TANGENTIAL_COMPLEX_CONFIG_H_

#include <cstddef>

// ========================= Debugging & profiling =============================
// #define GUDHI_TC_PROFILING
// #define GUDHI_TC_VERY_VERBOSE
// #define GUDHI_TC_PERFORM_EXTRA_CHECKS
// #define GUDHI_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES

// ========================= Strategy ==========================================
#define GUDHI_TC_PERTURB_POSITION
// #define GUDHI_TC_PERTURB_WEIGHT

// ========================= Parameters ========================================

// PCA will use GUDHI_TC_BASE_VALUE_FOR_PCA^intrinsic_dim points
const std::size_t GUDHI_TC_BASE_VALUE_FOR_PCA = 5;

#endif  // TANGENTIAL_COMPLEX_CONFIG_H_
