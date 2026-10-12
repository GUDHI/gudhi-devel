/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PYTHON_INTERFACES_CGAL_RANDOM_H_
#define PYTHON_INTERFACES_CGAL_RANDOM_H_

// Must be done before #include <CGAL/Random.h>
namespace CGAL {
  class Random;
#ifdef _WIN32
  #if defined GUDHI_DEFAULT_CGAL_RANDOM_DLL_IMPORT
    extern __declspec(dllimport) Random& get_default_random();
  #endif
  #if defined GUDHI_DEFAULT_CGAL_RANDOM_DLL_EXPORT
    __declspec(dllexport) Random& get_default_random();
  #endif
#else
  Random& get_default_random() __attribute__((visibility("default")));
#endif
}

#include <CGAL/Random.h>  // for CGAL::get_default_random()


#endif  // PYTHON_INTERFACES_CGAL_RANDOM_H_
