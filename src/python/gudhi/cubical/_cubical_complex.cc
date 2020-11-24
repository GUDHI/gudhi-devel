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


template<typename CubicalComplexOptions = Bitmap_cubical_complex_base<double>>
class Cubical_complex_interface_ : public Bitmap_cubical_complex<CubicalComplexOptions> {
 public:
  Cubical_complex_interface_(const std::vector<unsigned>& dimensions,
                            const std::vector<double>& top_dimensional_cells)
  : Bitmap_cubical_complex<CubicalComplexOptions>(dimensions, top_dimensional_cells) {
  }

  Cubical_complex_interface_(const std::vector<unsigned>& dimensions,
                            const std::vector<double>& top_dimensional_cells,
                            const std::vector<bool>& periodic_dimensions)
  : Bitmap_cubical_complex<CubicalComplexOptions>(dimensions, top_dimensional_cells, periodic_dimensions) {
  }

  Cubical_complex_interface_(const std::string& perseus_file)
  : Bitmap_cubical_complex<CubicalComplexOptions>(perseus_file.c_str()) {
  }

  std::size_t get_dimension() const { return Bitmap_cubical_complex<CubicalComplexOptions>::dimension(); }
};

using Persistent_cohomology_cubical_interface_ = Gudhi::Persistent_cohomology_interface<Cubical_complex_interface_<>>;

using Periodic_cubical_complex_interface_ =
    Cubical_complex_interface_<Gudhi::cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double>>;
using Persistent_cohomology_periodic_cubical_interface_ =
    Gudhi::Persistent_cohomology_interface<Periodic_cubical_complex_interface_>;

PYBIND11_MODULE(_cubical_complex, m) {
  // Cubical complex
  py::class_<Cubical_complex_interface_<>>(m, "Cubical_complex_interface_")
    .def(py::init<const std::string &>())
    .def(py::init<const std::vector<unsigned>&, const std::vector<double>&>())
    .def("num_simplices", &Cubical_complex_interface_<>::num_simplices)
    .def("dimension", &Cubical_complex_interface_<>::get_dimension);

  py::class_<Persistent_cohomology_cubical_interface_>(m, "Persistent_cohomology_cubical_interface_")
    .def(py::init<Cubical_complex_interface_<>*, bool>())
    .def("compute_persistence", &Persistent_cohomology_cubical_interface_::compute_persistence)
    .def("get_persistence", &Persistent_cohomology_cubical_interface_::get_persistence)
    .def("cofaces_of_cubical_persistence_pairs", &Persistent_cohomology_cubical_interface_::cofaces_of_cubical_persistence_pairs)
    .def("betti_numbers", &Persistent_cohomology_cubical_interface_::betti_numbers)
    .def("persistent_betti_numbers", &Persistent_cohomology_cubical_interface_::persistent_betti_numbers)
    .def("intervals_in_dimension", &Persistent_cohomology_cubical_interface_::intervals_in_dimension);

  // Periodic cubical complex
  py::class_<Periodic_cubical_complex_interface_>(m, "Periodic_cubical_complex_interface_")
    .def(py::init<const std::string &>())
    .def(py::init<const std::vector<unsigned>&, const std::vector<double>&, const std::vector<bool>&>())
    .def("num_simplices", &Periodic_cubical_complex_interface_::num_simplices)
    .def("dimension", &Periodic_cubical_complex_interface_::get_dimension);

  py::class_<Persistent_cohomology_periodic_cubical_interface_>(m, "Persistent_cohomology_periodic_cubical_interface_")
    .def(py::init<Periodic_cubical_complex_interface_*, bool>())
    .def("compute_persistence", &Persistent_cohomology_periodic_cubical_interface_::compute_persistence)
    .def("get_persistence", &Persistent_cohomology_periodic_cubical_interface_::get_persistence)
    .def("cofaces_of_cubical_persistence_pairs",
        &Persistent_cohomology_periodic_cubical_interface_::cofaces_of_cubical_persistence_pairs)
    .def("betti_numbers", &Persistent_cohomology_periodic_cubical_interface_::betti_numbers)
    .def("persistent_betti_numbers", &Persistent_cohomology_periodic_cubical_interface_::persistent_betti_numbers)
    .def("intervals_in_dimension", &Persistent_cohomology_periodic_cubical_interface_::intervals_in_dimension);
}
