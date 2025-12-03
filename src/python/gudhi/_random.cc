/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <nanobind/nanobind.h>

#include <python_interfaces/random_utils.h>

namespace nb = nanobind;

NB_MODULE(_random_ext, m) {
  // Based on https://doc.cgal.org/latest/Generator/index.html
  m.attr("__license__") = "LGPL v3";
  
  nb::class_<Gudhi::random::Random_generator> (m, "GudhiRandomGenerator")
      .def(nb::init<>())
      .def(nb::init<long>())
      .def("setup_bitgen", &Gudhi::random::setup_bitgen);
}
