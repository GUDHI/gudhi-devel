/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <nanobind/nanobind.h>

// For Windows, where _default_random is the dll provider for the other clients
// Must be done before #include <gudhi/random.h>
#define GUDHI_DEFAULT_RANDOM_DLL_IMPORT
#include <gudhi/random.h>

NB_MODULE(_default_random_seed_ext, m) {
  m.attr("__license__") = "MIT";
  m.def("_set_seed", [](int seed) { Gudhi::random::set_seed(seed); });
}
