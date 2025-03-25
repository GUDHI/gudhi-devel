/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 */

#include <../include/Simplex_tree_interface.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb=nanobind;

NB_MODULE(simplex_tree, m) {
      m.attr("__license__") = "GPL v3";
      nb::class_<Gudhi::Simplex_tree_interface>(m, "Simplex_tree")
          .def(nb::init<>())
          .def("find_simplex",
               &Gudhi::Simplex_tree_interface::find_simplex,
               nb::arg("simplex"),
               R"pbdoc()pbdoc")
          ;

}
