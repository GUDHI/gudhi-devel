/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2017 Inria
 *
 *    Modification(s):
 *      - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings
 *      - YYYY/MM Author: Description of the modification
 */

// /////////////////////////////////////////////////////////////////////////////
// Reader_utils_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/map.h>

#include <gudhi/reader_utils.h>

namespace nb = nanobind;

NB_MODULE(_reader_utils_ext, m)
{
  m.attr("__license__") = "MIT";

  m.def("read_matrix_from_csv_file",
        &Gudhi::read_lower_triangular_matrix_from_csv_file<double>,
        nb::arg("filename"),
        nb::arg("separator") = ';',
        "read_matrix_from_csv_file");
  m.def("read_pers_intervals_grouped_by_dimension",
        &Gudhi::read_persistence_intervals_grouped_by_dimension,
        nb::arg("filename"),
        "read_pers_intervals_grouped_by_dimension");
  m.def("read_pers_intervals_in_dimension",
        &Gudhi::read_persistence_intervals_in_dimension,
        nb::arg("filename"),
        nb::arg("only_this_dim") = -1,
        "read_pers_intervals_in_dimension");
}

//
// _reader_utils.cc ends here