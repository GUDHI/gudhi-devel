/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
*    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thibaud Kloczko
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include "Reader_utils_interface.h"

// /////////////////////////////////////////////////////////////////////////////
// Reader_utils_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/map.h>

#include "../../Bitmap_cubical_complex/doc/Gudhi_Cubical_Complex_doc.h"

namespace nb = nanobind;

NB_MODULE(_reader_utils_ext, m) {
 m.attr("__license__") = "GPL v3";

 m.def("read_matrix_from_csv_file", &Gudhi::read_matrix_from_csv_file, "read_matrix_from_csv_file");
 m.def("read_pers_intervals_grouped_by_dimension", &Gudhi::read_pers_intervals_grouped_by_dimension, "read_pers_intervals_grouped_by_dimension");
 m.def("read_pers_intervals_in_dimension", &Gudhi::read_pers_intervals_in_dimension, "read_pers_intervals_in_dimension");
}

//
// _reader_utils.cc ends here