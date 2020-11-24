/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>

#include <Persistent_cohomology_interface.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <iostream>
#include <vector>
#include <string>
#include <memory>       // for std::unique_ptr

namespace py = pybind11;
using namespace Gudhi::cubical_complex;

using Base = Gudhi::cubical_complex::Bitmap_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>>;

class Cubical_complex_interface_ : public Base {
 public:
  Cubical_complex_interface_(const std::vector<unsigned>& dimensions,
                            const std::vector<double>& top_dimensional_cells)
  : Base(dimensions, top_dimensional_cells) {
  }

  Cubical_complex_interface_(const std::vector<unsigned>& dimensions,
                            const std::vector<double>& top_dimensional_cells,
                            const std::vector<bool>& periodic_dimensions)
  : Base(dimensions, top_dimensional_cells, periodic_dimensions) {
  }

  Cubical_complex_interface_(const std::string& perseus_file)
  : Base(perseus_file.c_str()) {
  }

  std::size_t get_dimension() const { return dimension(); }

  //Gudhi::Persistent_cohomology_interface<Cubical_complex_interface_> 
  void get_persistence() {
    Gudhi::Persistent_cohomology_interface<Cubical_complex_interface_> persistence(this, true);
    //return persistence;
  }
};

using Persistent_cohomology_cubical_interface_ = Gudhi::Persistent_cohomology_interface<Cubical_complex_interface_>;

PYBIND11_MODULE(_cubical_complex, m) {
  py::class_<Cubical_complex_interface_>(m, "Cubical_complex_interface_")
    .def(py::init<const std::string &>())
    .def(py::init<const std::vector<unsigned>&, const std::vector<double>&>())
    .def("num_simplices", &Cubical_complex_interface_::num_simplices)
    .def("dimension", &Cubical_complex_interface_::get_dimension)
    .def("get_persistence", &Cubical_complex_interface_::get_persistence);

  py::class_<Persistent_cohomology_cubical_interface_>(m, "Persistent_cohomology_cubical_interface_")
    .def(py::init<Cubical_complex_interface_*, bool>())
    .def("compute_persistence", &Persistent_cohomology_cubical_interface_::compute_persistence)
    .def("get_persistence", &Persistent_cohomology_cubical_interface_::get_persistence)
    .def("cofaces_of_cubical_persistence_pairs", &Persistent_cohomology_cubical_interface_::cofaces_of_cubical_persistence_pairs)
    .def("betti_numbers", &Persistent_cohomology_cubical_interface_::betti_numbers)
    .def("persistent_betti_numbers", &Persistent_cohomology_cubical_interface_::persistent_betti_numbers)
    .def("intervals_in_dimension", &Persistent_cohomology_cubical_interface_::intervals_in_dimension);
}