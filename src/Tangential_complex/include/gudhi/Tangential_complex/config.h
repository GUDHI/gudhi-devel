/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
