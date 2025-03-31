/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
*    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Jean Luc Szpyrka
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include "Subsampling_interface.h"

// /////////////////////////////////////////////////////////////////////////////
// Reader_utils_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

NB_MODULE(_subsampling_ext, m) {
 m.attr("__license__") = "MIT (GPL v3 for sparsify_point_set)";

 m.def("subsampling_n_farthest_points",&Gudhi::subsampling::subsampling_n_farthest_points_from_file, "");
 m.def("subsampling_n_farthest_points_from_file", &Gudhi::subsampling::subsampling_n_farthest_points_from_file, "");
 m.def("subsampling_n_random_points", &Gudhi::subsampling::subsampling_n_random_points, "");
 m.def("subsampling_n_random_points_from_file", &Gudhi::subsampling::subsampling_n_random_points_from_file, "");
 m.def("subsampling_sparsify_points", &Gudhi::subsampling::subsampling_sparsify_points, "");
 m.def("subsampling_sparsify_points_from_file", &Gudhi::subsampling::subsampling_sparsify_points_from_file, "");
}

//
// _subsampling.cc ends here
