/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

// For Windows, where _default_random is the dll provider for the other clients
// Must be done before #include <gudhi/random.h>
#define GUDHI_DEFAULT_RANDOM_DLL_EXPORT
#include <gudhi/random.h>
